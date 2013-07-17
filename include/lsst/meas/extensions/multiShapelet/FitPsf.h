// -*- lsst-c++ -*-
/* 
 * LSST Data Management System
 * Copyright 2012 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
#ifndef MULTISHAPELET_FitPsf_h_INCLUDED
#define MULTISHAPELET_FitPsf_h_INCLUDED

#include "ndarray.h"

#include "lsst/shapelet.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/Algorithm.h"
#include "lsst/meas/extensions/multiShapelet/MultiGaussian.h"
#include "lsst/meas/extensions/multiShapelet/ModelInputHandler.h"
#include "lsst/meas/extensions/multiShapelet/HybridOptimizer.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

class FitPsfAlgorithm;

class MultiGaussianObjective;

class FitPsfControl : public algorithms::AlgorithmControl {
public:

    LSST_CONTROL_FIELD(innerOrder, int, "Shapelet order of inner expansion (0 == Gaussian)");
    LSST_CONTROL_FIELD(outerOrder, int, "Shapelet order of outer expansion (0 == Gaussian)");
    LSST_CONTROL_FIELD(minRadius, double, "Minimum inner radius in pixels.");
    LSST_CONTROL_FIELD(minAxisRatio, double, "Minimum axis ratio for ellipse (b/a).");
    LSST_CONTROL_FIELD(radiusRatio, double, "outer radius divided by inner radius (fixed)");
    LSST_CONTROL_FIELD(peakRatio, double, 
                       "outer Gaussian peak height divided by inner Gaussian peak height;"
                       " held fixed in double-Gaussian ellipse fit, then allowed to vary"
                       " when shapelets coefficients are fit and ellipses are held fixed."
    );
    LSST_CONTROL_FIELD(initialRadius, double, "Initial radius of inner component in pixels");
    LSST_CONTROL_FIELD(useApproximateExp, bool, "Use fast approximate exponential (good to ~1E-4)");

    PTR(FitPsfControl) clone() const { return boost::static_pointer_cast<FitPsfControl>(_clone()); }

    PTR(FitPsfAlgorithm) makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata = PTR(daf::base::PropertyList)(),
        algorithms::AlgorithmMap const & others = algorithms::AlgorithmMap(),
        bool isForced = false
    ) const;
    
    MultiGaussian getMultiGaussian() const;

    FitPsfControl() :
        algorithms::AlgorithmControl("multishapelet.psf", 2.0),
        innerOrder(2), outerOrder(2), minRadius(0.1), minAxisRatio(0.1),
        radiusRatio(2.0), peakRatio(0.1), initialRadius(1.5), useApproximateExp(false)
    {}

private:

    using algorithms::AlgorithmControl::_makeAlgorithm;

    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

/**
 *  @brief A Multi-Shapelet model for a local PSF.
 *
 *  This is an extension of an elliptical double-Gaussian, where each Gaussian is replaced
 *  with a Gauss-Hermite expansion of arbitrary order (when both orders are zero, it is
 *  exactly a double-Gaussian).  The inner and outer components have the same basis
 *  ellipticity, but the ratio of their radii is not fixed.
 *
 *  See lsst::shapelet for the precise definitions of the basis functions.
 */
struct FitPsfModel {

    ndarray::Array<double,1,1> inner; ///< shapelet coefficients of inner expansion
    ndarray::Array<double,1,1> outer; ///< shapelet coefficients of outer expansion
    afw::geom::ellipses::Quadrupole ellipse; ///< ellipse corresponding to inner expansion
    double radiusRatio; ///< radius of outer expansion divided by radius of inner expansion (fixed)
    double chisq; ///< reduced chi^2
    bool failedMaxIter; ///< set to true if the optimizer hit the maximum number of iterations
    bool failedTinyStep; ///< set to true if the optimizer step size got too small to make progress
    bool failedMinRadius; ///< set to true if the best-fit radius was at the minimum constraint
    bool failedMinAxisRatio; ///< set to true if the best-fit axis ratio was at the minimum constraint

    bool hasFailed() const { return failedMaxIter || failedTinyStep || failedMinRadius || failedMinAxisRatio; }

    /**
     *  @brief Construct a model from a double-Gaussian optimization parameter vector.
     *
     *  This is used internally by FitPsfAlgorithm to construct the model after fitting
     *  an elliptical double-Gaussian, and the 3-element parameter vector should correspond
     *  to the parameters used by FitPsfAlgorithm::makeObjective and
     *  FitPsfAlgorithm::makeOptimizer.  By using those functions and this constructor,
     *  we can inspect the optimizer at each step in Python.
     */
    FitPsfModel(
        FitPsfControl const & ctrl,
        double amplitude,
        ndarray::Array<double const,1,1> const & parameters
    );

    /// @brief Construct by extracting saved values from a Record.
    FitPsfModel(FitPsfControl const & ctrl, afw::table::BaseRecord const & source);

    /// @brief Deep copy constructor.
    FitPsfModel(FitPsfModel const & other);

    /// @brief Deep assignment operator.
    FitPsfModel & operator=(FitPsfModel const & other);
    
    MultiGaussian getMultiGaussian() const;

    /**
     *  @brief Return a MultiShapeletFunction representation of the model.
     *
     *  The elements will be [inner, outer].
     */
    shapelet::MultiShapeletFunction asMultiShapelet(
        afw::geom::Point2D const & center = afw::geom::Point2D()
    ) const;

};

class FitPsfAlgorithm : public algorithms::Algorithm {
public:

    typedef FitPsfControl Control;
    typedef FitPsfModel Model;

    /// @brief Construct an algorithm instance and register its fields with a Schema.
    FitPsfAlgorithm(FitPsfControl const & ctrl, afw::table::Schema & schema);

    /// @brief Return the control object
    FitPsfControl const & getControl() const {
        return static_cast<FitPsfControl const &>(algorithms::Algorithm::getControl());
    }

    /**
     *  @brief Return an Objective that can be used to fit an elliptical double-Gaussian to the image.
     *
     *  This is provided primarily for testing and debugging purposes.
     */
    static PTR(MultiGaussianObjective) makeObjective(
        FitPsfControl const & ctrl,
        ModelInputHandler const & inputs
    );

    /**
     *  @brief Return an optimizer that can be used to fit an elliptical double-Gaussian to the image.
     *
     *  This is provided primarily for testing and debugging purposes; the user can create an optimizer,
     *  step through it, and use the FitPsfModel constructor that takes a parameter vector to
     *  visualize its progress.
     */
    static HybridOptimizer makeOptimizer(
        FitPsfControl const & ctrl,
        ModelInputHandler const & inputs
    );

    /**
     *  @brief Given a double-Gaussian-only model, fit additional shapelet terms.
     *
     *  This is provided primarily for testing and debugging purposes; it is the second
     *  part of apply(), after the nonlinear double-Gaussian fit.
     */
    static void fitShapeletTerms(
        FitPsfControl const & ctrl,
        ModelInputHandler const & inputs,
        FitPsfModel & model
    );

    /**
     *  @brief Fit to a local PSF or kernel image.
     *
     *  @param[in] ctrl           Details of the model to fit.
     *  @param[in] inputs         Inputs that determine the data to be fit.
     */
    static FitPsfModel apply(
        FitPsfControl const & ctrl,
        ModelInputHandler const & inputs
    );

    /**
     *  @brief Fit a PSF object evaluated at a point.
     *
     *  @param[in] ctrl           Details of the model to fit.
     *  @param[in] psf            PSF object
     *  @param[in] center         Point at which to evaluate the PSF.
     */
    static FitPsfModel apply(
        FitPsfControl const & ctrl,
        afw::detection::Psf const & psf,
        afw::geom::Point2D const & center
    );

    /**
     *  @brief Fit a PSF object evaluated at a point, returning a FitPsfModel
     *         and saving the to the given record.
     *
     *  @param[in] record         Record in which to save results
     *  @param[in] psf            PSF object
     *  @param[in] center         Point at which to evaluate the PSF
     */
    FitPsfModel apply(
        afw::table::BaseRecord & record,
        afw::detection::Psf const & psf,
        afw::geom::Point2D const & center
    ) const;

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(FitPsfAlgorithm);

    afw::table::Key< afw::table::Array<float> > _innerKey;
    afw::table::Key< afw::table::Array<float> > _outerKey;
    afw::table::Key< afw::table::Moments<float> > _ellipseKey;
    afw::table::Key<float> _chisqKey;
    afw::table::Key<float> _integralKey;
    afw::table::Key< afw::table::Flag > _flagKey;
    afw::table::Key< afw::table::Flag > _flagMaxIterKey;
    afw::table::Key< afw::table::Flag > _flagTinyStepKey;
    afw::table::Key< afw::table::Flag > _flagMinRadiusKey;
    afw::table::Key< afw::table::Flag > _flagMinAxisRatioKey;
};

inline PTR(FitPsfAlgorithm) FitPsfControl::makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata,
    algorithms::AlgorithmMap const & others,
    bool isForced
) const {
    return boost::static_pointer_cast<FitPsfAlgorithm>(_makeAlgorithm(schema, metadata, others, isForced));
}

}}}} // namespace lsst::meas::extensions::multiShapelet

#endif // !MULTISHAPELET_FitPsf_h_INCLUDED
