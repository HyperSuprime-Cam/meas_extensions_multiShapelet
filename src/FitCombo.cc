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

#include "lsst/utils/ieee.h"
#include "lsst/meas/extensions/multiShapelet/FitCombo.h"
#include "lsst/shapelet/ModelBuilder.h"
#include "lsst/afw/math/LeastSquares.h"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"


namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

//------------ FitComboControl ----------------------------------------------------------------------------

PTR(algorithms::AlgorithmControl) FitComboControl::_clone() const {
    return boost::make_shared<FitComboControl>(*this);
}

PTR(algorithms::Algorithm) FitComboControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metxadata,
    algorithms::AlgorithmMap const & others,
    bool isForced
) const {
    return boost::make_shared<FitComboAlgorithm>(*this, boost::ref(schema), others);
}

//------------ FitComboModel ------------------------------------------------------------------------------

FitComboModel::FitComboModel(FitComboControl const & ctrl) :
    devFrac(0.0),
    flux(std::numeric_limits<double>::quiet_NaN()), fluxErr(std::numeric_limits<double>::quiet_NaN())
{}

FitComboModel::FitComboModel(FitComboModel const & other) :
    devFrac(other.devFrac), flux(other.flux), fluxErr(other.fluxErr)
{}

FitComboModel & FitComboModel::operator=(FitComboModel const & other) {
    if (&other != this) {
        devFrac = other.devFrac;
        flux = other.flux;
        fluxErr = other.fluxErr;
    }
    return *this;
}

//------------ FitComboAlgorithm --------------------------------------------------------------------------

namespace {

// helper function for constructor to get control objects for dependencies
template <typename Ctrl>
CONST_PTR(Ctrl) getDependency(algorithms::AlgorithmMap const & others, std::string const & name) {
    algorithms::AlgorithmMap::const_iterator i = others.find(name);
    if (i == others.end()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            (boost::format("Algorithm with name '%s' not found; needed by FitCombo.") % name).str()
        );
    }
    CONST_PTR(Ctrl) result = boost::dynamic_pointer_cast<Ctrl const>(i->second->getControl().clone());
    if (!result) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            (boost::format("Algorithm with name '%s' does not have the correct type.") % name).str()
        );
    }
    return result;
}

} // anonymous

FitComboAlgorithm::FitComboAlgorithm(
    FitComboControl const & ctrl,
    afw::table::Schema & schema,
    algorithms::AlgorithmMap const & others
) :
    algorithms::Algorithm(ctrl),
    _fluxKeys(
        afw::table::addFluxFields(
            schema, ctrl.name + ".flux",
            "combined flux of the linear combination of fixed-profile models"
        )
    ),
    _devFracKey(
        schema.addField<float>(
            ctrl.name + ".devfrac",
            "fraction of total flux in the de Vaucouleur component"
        )),
    _expComponentCtrl(getDependency<FitProfileControl>(others, ctrl.expName)),
    _devComponentCtrl(getDependency<FitProfileControl>(others, ctrl.devName)),
    _psfCtrl(getDependency<FitPsfControl>(others, ctrl.psfName))
{}

template <typename PixelT>
ModelInputHandler FitComboAlgorithm::adjustInputs(
    FitComboControl const & ctrl,
    FitPsfModel const & psfModel,
    FitProfileModel const & expComponent,
    FitProfileModel const & devComponent,
    afw::detection::Footprint const & footprint,
    afw::image::MaskedImage<PixelT> const & image,
    afw::geom::Point2D const & center
) {
    afw::image::MaskPixel badPixelMask(0);
        badPixelMask = afw::image::Mask<>::getPlaneBitMask(ctrl.badMaskPlanes);
    if (ctrl.radiusInputFactor > 0.0) {
        std::vector<afw::geom::ellipses::Ellipse> boundsEllipses;
        boundsEllipses.push_back(afw::geom::ellipses::Ellipse(expComponent.ellipse, center));
        boundsEllipses.back().getCore().scale(ctrl.radiusInputFactor);
        boundsEllipses.push_back(afw::geom::ellipses::Ellipse(devComponent.ellipse, center));
        boundsEllipses.back().getCore().scale(ctrl.radiusInputFactor);
        return ModelInputHandler(image, center,
                                 boundsEllipses, footprint, ctrl.growFootprint, 
                                 badPixelMask, ctrl.usePixelWeights);
    } else {
        return ModelInputHandler(image, center, footprint, ctrl.growFootprint, 
                                 badPixelMask, ctrl.usePixelWeights);
    }
}

namespace {

// helper function for FitComboAlgorithm::apply
void buildComponentModel(
    FitProfileModel const & component,
    FitPsfModel const & psfModel,
    shapelet::ModelBuilder<double> & builder,
    ModelInputHandler const & inputs,
    ndarray::Array<double,1,1> const & output,
    double weight = 1.0
) {
    typedef shapelet::MultiShapeletFunction MSF;
   MSF msf = component.asMultiShapelet(afw::geom::Point2D())
       .convolve(psfModel.asMultiShapelet());
   msf.normalize(weight);
   for (
       MSF::ElementList::const_iterator i = msf.getElements().begin();
       i != msf.getElements().end();
       ++i
   ) {
       builder.update(i->getEllipse().getCore());
       builder.addModelVector(i->getOrder(), i->getCoefficients(), output);
   }
}

} // anonymous

FitComboModel FitComboAlgorithm::apply(
    FitComboControl const & ctrl,
    FitPsfModel const & psfModel,
    FitProfileModel const & expComponent,
    FitProfileModel const & devComponent,
    ModelInputHandler const & inputs
) {
    FitComboModel model(ctrl);
    typedef shapelet::MultiShapeletFunction MSF;
    ndarray::Array<double,2,2> matrixT = ndarray::allocate(2, inputs.getSize());
    shapelet::ModelBuilder<double> builder(inputs.getX(), inputs.getY(), ctrl.useApproximateExp);
    ndarray::Array<double,2,-2> matrix(matrixT.transpose());
    matrixT.deep() = 0.0;
    buildComponentModel(expComponent, psfModel, builder, inputs, matrixT[0]);
    if (!inputs.getWeights().isEmpty()) {
        matrixT[0].asEigen<Eigen::ArrayXpr>() *= inputs.getWeights().asEigen<Eigen::ArrayXpr>();
    }
    buildComponentModel(devComponent, psfModel, builder, inputs, matrixT[1]);
    if (!inputs.getWeights().isEmpty()) {
        matrixT[1].asEigen<Eigen::ArrayXpr>() *= inputs.getWeights().asEigen<Eigen::ArrayXpr>();
    }
    // We should really do constrained linear least squares to get the errors right, but this
    // produces the same result for the fluxes, and we don't have a constrained solver handy.
    afw::math::LeastSquares lstsq = afw::math::LeastSquares::fromDesignMatrix(matrix, inputs.getData());
    if (lstsq.getSolution()[0] < 0.0) {
        if (lstsq.getSolution()[1] < 0.0) {
            throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "measured negative flux");
        }
        model.devFrac = 1.0;
        model.flux = devComponent.flux;
        model.fluxErr = devComponent.fluxErr;
    } else if (lstsq.getSolution()[1] < 0.0) {
        model.devFrac = 0.0;
        model.flux = expComponent.flux;
        model.fluxErr = expComponent.fluxErr;
    } else {
        model.flux = lstsq.getSolution().asEigen().sum();
        model.devFrac = lstsq.getSolution()[1] / model.flux;
#if 0 // don't know why this doesn't work; numbers are way too big
        model.fluxErr = std::sqrt(
            lstsq.getSolution().asEigen().dot(
                lstsq.getCovariance().asEigen() * lstsq.getSolution().asEigen()
            )
        );
#else // this is incorrect, but a good-enough workaround for now: weighted average in quadrature
        model.fluxErr = std::sqrt(expComponent.fluxErr * expComponent.fluxErr * (1.0 - model.devFrac)
                                  + devComponent.fluxErr * devComponent.fluxErr * model.devFrac);
#endif
    }
    return model;
}

template <typename PixelT>
void FitComboAlgorithm::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    source.set(_fluxKeys.flag, true);
    if (!exposure.hasPsf()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Cannot run FitComboAlgorithm without a PSF."
        );
    }
    FitPsfModel psfModel(*_psfCtrl, source);
    FitProfileModel expComponent(*_expComponentCtrl, source);
    FitProfileModel devComponent(*_devComponentCtrl, source);
    if (expComponent.fluxFlag || devComponent.fluxFlag) {
        return; // Don't bother trying linear components if one of the inputs failed.
    }
    assert(lsst::utils::isfinite(expComponent.ellipse.getArea()));
    assert(lsst::utils::isfinite(devComponent.ellipse.getArea()));
    ModelInputHandler inputs = adjustInputs(
        getControl(), psfModel, expComponent, devComponent, *source.getFootprint(),
        exposure.getMaskedImage(), center
    );
    FitComboModel model = apply(getControl(), psfModel, expComponent, devComponent, inputs);

    source.set(_devFracKey, model.devFrac);
    source.set(_fluxKeys.meas, model.flux);
    source.set(_fluxKeys.err, model.fluxErr);
    source.set(_fluxKeys.flag, false);
}


template <typename PixelT>
void FitComboAlgorithm::_applyForced(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center,
    afw::table::SourceRecord const & reference,
    afw::geom::AffineTransform const & refToMeas
) const {
    source.set(_fluxKeys.flag, true);
    if (!exposure.hasPsf()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Cannot run FitComboAlgorithm without a PSF."
        );
    }
    FitPsfModel psfModel(*_psfCtrl, source);
    FitProfileModel expComponent(*_expComponentCtrl, source);
    FitProfileModel devComponent(*_devComponentCtrl, source);
    if (expComponent.fluxFlag || devComponent.fluxFlag) {
        return; // Don't bother trying linear components if one of the inputs failed.
    }
    assert(lsst::utils::isfinite(expComponent.ellipse.getArea()));
    assert(lsst::utils::isfinite(devComponent.ellipse.getArea()));
    afw::table::SubSchema s = reference.getSchema()[getControl().name];
    ModelInputHandler inputs = adjustInputs(
        getControl(), psfModel, expComponent, devComponent, *source.getFootprint(),
        exposure.getMaskedImage(), center
    );
    FitComboModel model(getControl());
    model.devFrac = reference.get(s.find<float>("devfrac").key);
    shapelet::ModelBuilder<double> builder(inputs.getX(), inputs.getY());
    ndarray::Array<double,1,1> matrix = ndarray::allocate(inputs.getSize());
    matrix.deep() = 0.0;
    buildComponentModel(expComponent, psfModel, builder, inputs, matrix, 1.0 - model.devFrac);
    buildComponentModel(devComponent, psfModel, builder, inputs, matrix, model.devFrac);
    if (!inputs.getWeights().isEmpty()) {
        matrix.asEigen<Eigen::ArrayXpr>() *= inputs.getWeights().asEigen<Eigen::ArrayXpr>();
    }
    // 1-d linear least squares...
    double a = matrix.asEigen<Eigen::MatrixXpr>().squaredNorm();
    double b = matrix.asEigen<Eigen::MatrixXpr>().dot(inputs.getData().asEigen<Eigen::MatrixXpr>());
    model.flux = b / a;
    model.fluxErr = std::sqrt(1.0 / a);

    source.set(_devFracKey, model.devFrac);
    source.set(_fluxKeys.meas, model.flux);
    source.set(_fluxKeys.err, model.fluxErr);
    source.set(_fluxKeys.flag, false);
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(FitComboAlgorithm);

template
ModelInputHandler FitComboAlgorithm::adjustInputs(
    FitComboControl const & ctrl,
    FitPsfModel const & psfModel,
    FitProfileModel const & expComponent,
    FitProfileModel const & devComponent,
    afw::detection::Footprint const & footprint,
    afw::image::MaskedImage< float > const & image,
    afw::geom::Point2D const & center
						  );

template
ModelInputHandler FitComboAlgorithm::adjustInputs(
    FitComboControl const & ctrl,
    FitPsfModel const & psfModel,
    FitProfileModel const & expComponent,
    FitProfileModel const & devComponent,
    afw::detection::Footprint const & footprint,
    afw::image::MaskedImage< double > const & image,
    afw::geom::Point2D const & center
						  );

}}}} // namespace lsst::meas::extensions::multiShapelet
