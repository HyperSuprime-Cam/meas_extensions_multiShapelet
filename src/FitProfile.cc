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
#include "lsst/meas/extensions/multiShapelet/FitProfile.h"
#include "lsst/meas/extensions/multiShapelet/MultiGaussianObjective.h"
#include "lsst/meas/extensions/multiShapelet/MultiGaussianRegistry.h"
#include "lsst/shapelet/ModelBuilder.h"
#include "lsst/afw/math/LeastSquares.h"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

//------------ FitProfileControl ----------------------------------------------------------------------------

PTR(algorithms::AlgorithmControl) FitProfileControl::_clone() const {
    return boost::make_shared<FitProfileControl>(*this);
}

PTR(algorithms::Algorithm) FitProfileControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata,
    algorithms::AlgorithmMap const & others,
    bool isForced
) const {
    return boost::make_shared<FitProfileAlgorithm>(*this, boost::ref(schema), others, isForced);
}

//------------ FitProfileModel ------------------------------------------------------------------------------

FitProfileModel::FitProfileModel(
    FitProfileControl const & ctrl,
    double amplitude,
    ndarray::Array<double const,1,1> const & parameters
) :
    profile(ctrl.profile), flux(amplitude), fluxErr(0.0),
    ellipse(MultiGaussianObjective::EllipseCore(parameters[0], parameters[1], parameters[2])),
    chisq(std::numeric_limits<double>::quiet_NaN()), fluxFlag(false),
    flagMaxIter(false), flagTinyStep(false), flagMinRadius(false), flagMinAxisRatio(false),
    flagLargeArea(false)
{}

FitProfileModel::FitProfileModel(
    FitProfileControl const & ctrl, afw::table::SourceRecord const & source,
    bool loadPsfFactorModel
) :
    profile(ctrl.profile), flux(1.0), fluxErr(0.0), ellipse(),
    chisq(std::numeric_limits<double>::quiet_NaN()), fluxFlag(false),
    flagMaxIter(false), flagTinyStep(false), flagMinRadius(false), flagMinAxisRatio(false),
    flagLargeArea(false)
{
    afw::table::SubSchema s = source.getSchema()[ctrl.name];
    if (loadPsfFactorModel) {
        flux = source.get(s.find< float >("psffactor").key);
        fluxFlag = source.get(s.find<afw::table::Flag>("flags.psffactor").key);
        ellipse = source.get(s.find< afw::table::Moments<double> >("psffactor.ellipse").key);
    } else {
        flux = source.get(s.find< double >("flux").key);
        fluxErr = source.get(s.find< double >("flux.err").key);
        fluxFlag = source.get(s.find<afw::table::Flag>("flux.flags").key);
        ellipse = source.get(s.find< afw::table::Moments<double> >("ellipse").key);
        chisq = source.get(s.find<float>("chisq").key);
        flagMaxIter = source.get(s.find<afw::table::Flag>("flags.maxiter").key);
        flagTinyStep = source.get(s.find<afw::table::Flag>("flags.tinystep").key);
        flagLargeArea = source.get(s.find<afw::table::Flag>("flags.largearea").key);
    }
    assert(fluxFlag || lsst::utils::isfinite(ellipse.getArea()));
}

FitProfileModel::FitProfileModel(FitProfileModel const & other) :
    profile(other.profile), flux(other.flux), fluxErr(other.fluxErr), ellipse(other.ellipse),
    chisq(other.chisq),
    fluxFlag(other.fluxFlag),
    flagMaxIter(other.flagMaxIter),
    flagTinyStep(other.flagTinyStep),
    flagMinRadius(other.flagMinRadius),
    flagMinAxisRatio(other.flagMinAxisRatio),
    flagLargeArea(other.flagLargeArea)
{}

FitProfileModel & FitProfileModel::operator=(FitProfileModel const & other) {
    if (&other != this) {
        profile = other.profile;
        flux = other.flux;
        fluxErr = other.fluxErr;
        fluxFlag = other.fluxFlag;
        ellipse = other.ellipse;
        chisq = other.chisq;
        flagMaxIter = other.flagMaxIter;
        flagTinyStep = other.flagTinyStep;
        flagMinRadius = other.flagMinRadius;
        flagMinAxisRatio = other.flagMinAxisRatio;
        flagLargeArea = other.flagLargeArea;
    }
    return *this;
}

shapelet::MultiShapeletFunction FitProfileModel::asMultiShapelet(
    afw::geom::Point2D const & center
) const {
    shapelet::MultiShapeletFunction::ElementList elements;
    MultiGaussian const & multiGaussian = MultiGaussianRegistry::lookup(profile);
    for (MultiGaussian::const_iterator i = multiGaussian.begin(); i != multiGaussian.end(); ++i) {
        afw::geom::ellipses::Ellipse fullEllipse(ellipse, center);
        elements.push_back(i->makeShapelet(fullEllipse));
        elements.back().getCoefficients().asEigen() *= flux;
    }
    return shapelet::MultiShapeletFunction(elements);
}

//------------ FitProfileAlgorithm --------------------------------------------------------------------------

FitProfileAlgorithm::FitProfileAlgorithm(
    FitProfileControl const & ctrl,
    afw::table::Schema & schema,
    algorithms::AlgorithmMap const & others,
    bool isForced
) :
    algorithms::Algorithm(ctrl),
    _fluxKeys(
        afw::table::addFluxFields(
            schema, ctrl.name + ".flux",
            "flux of multi-Gaussian approximation to the profile, integrated to infinite radius"
        )
    ),
    _ellipseKey(
        schema.addField< afw::table::Moments<double> >(
            ctrl.name + ".ellipse",
            "half-light radius ellipse"
        )),
    _psfEllipseKey(
        schema.addField< afw::table::Moments<double> >(
            ctrl.name + ".psffactor.ellipse",
            "half-light radius ellipse for fit to PSF model realization"
        )),
    _chisqKey(
        schema.addField<float>(
            ctrl.name + ".chisq",
            "reduced chi^2"
        )),
    _psfCtrl()
{
    if (!isForced) {
        _flagMaxIterKey = schema.addField<afw::table::Flag>(
            ctrl.name + ".flags.maxiter",
            "set if the optimizer ran into the maximum number of iterations limit"
        );
        _flagTinyStepKey = schema.addField<afw::table::Flag>(
            ctrl.name + ".flags.tinystep",
            "set if the optimizer step or trust region got so small no progress could be made"
        );
        _flagMinRadiusKey = schema.addField<afw::table::Flag>(
            ctrl.name + ".flags.constraint.r",
            "set if the best-fit radius was the minimum allowed by the constraint (not a failure)"
        );
        _flagMinAxisRatioKey = schema.addField<afw::table::Flag>(
            ctrl.name + ".flags.constraint.q",
            "set if the best-fit ellipticity was the maximum allowed by the constraint"
        );
        _flagLargeAreaKey = schema.addField<afw::table::Flag>(
            ctrl.name + ".flags.largearea",
            "set if the best-fit half-light ellipse area is larger than the number of pixels used"
        );
    }
    algorithms::AlgorithmMap::const_iterator i = others.find(ctrl.psfName);
    if (i == others.end()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            (boost::format("FitPsf with name '%s' not found; needed by FitProfile.") % ctrl.psfName).str()
        );
    }
    _psfCtrl = boost::dynamic_pointer_cast<FitPsfControl const>(i->second->getControl().clone());
    if (!_psfCtrl) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            (boost::format("Algorithm with name '%s' is not FitPsf.") % ctrl.psfName).str()
        );
    }
}

PTR(MultiGaussianObjective) FitProfileAlgorithm::makeObjective(
    FitProfileControl const & ctrl,
    FitPsfModel const & psfModel,
    ModelInputHandler const & inputs
) {
    return boost::make_shared<MultiGaussianObjective>(
        inputs, ctrl.getMultiGaussian(), psfModel.getMultiGaussian(), psfModel.ellipse,
        ctrl.minRadius, ctrl.minAxisRatio, ctrl.useApproximateExp
    );
}

HybridOptimizer FitProfileAlgorithm::makeOptimizer(
    FitProfileControl const & ctrl,
    FitPsfModel const & psfModel,
    MultiGaussianObjective::EllipseCore const & ellipse,
    ModelInputHandler const & inputs
) {
    PTR(Objective) obj = makeObjective(ctrl, psfModel, inputs);
    ndarray::Array<double,1,1> initial = ndarray::allocate(obj->getParameterSize());
    ellipse.writeParameters(initial.getData());
    HybridOptimizerControl optCtrl; // TODO: nest this in FitProfileControl
    optCtrl.tau = 1E-2;
    optCtrl.useCholesky = true;
    optCtrl.gTol = 1E-4;
    return HybridOptimizer(obj, initial, optCtrl);
}

template <typename PixelT>
ModelInputHandler FitProfileAlgorithm::adjustInputs(
    FitProfileControl const & ctrl,
    FitPsfModel const & psfModel,
    afw::geom::ellipses::Quadrupole & shape,
    afw::detection::Footprint const & footprint,
    afw::image::MaskedImage<PixelT> const & image,
    afw::geom::Point2D const & center,
    bool fixShape
) {
    MultiGaussianObjective::EllipseCore ellipse(shape);
    if (!(ellipse.getArea() > 0.0)) {  // phrasing comparison this way also guards against NaN
        ellipse = psfModel.ellipse;
        ellipse.scale(ctrl.minInitialRadius);
    } else {
        if (!fixShape && ctrl.deconvolveShape) {
            try {
                ellipse = ctrl.getMultiGaussian().deconvolve(
                    shape, psfModel.ellipse, psfModel.getMultiGaussian()
                );
            } catch (pex::exceptions::InvalidParameterException &) {
                ellipse = psfModel.ellipse;
                ellipse.scale(ctrl.minInitialRadius);
            }
        }
    }
    if (!fixShape) {
        // We never want to start with an ellipse smaller than the PSF or an ellipticity
        // on the constraint, because we might never find our way out.
        std::pair<bool,bool> constrained = MultiGaussianObjective::constrainEllipse(
            ellipse, psfModel.ellipse.getTraceRadius() * ctrl.minInitialRadius, ctrl.minAxisRatio
        );
        if (constrained.first || constrained.second) {
            ellipse = psfModel.ellipse;
            ellipse.scale(ctrl.minInitialRadius);
        }
        shape = ellipse;
    }
    afw::image::MaskPixel badPixelMask = afw::image::Mask<>::getPlaneBitMask(ctrl.badMaskPlanes);
    if (ctrl.radiusInputFactor > 0.0) {
        std::vector<afw::geom::ellipses::Ellipse> boundsEllipses;
        boundsEllipses.push_back(afw::geom::ellipses::Ellipse(shape, center));
        boundsEllipses.back().getCore().scale(ctrl.radiusInputFactor);
        return ModelInputHandler(image, center,
                                 boundsEllipses, footprint, ctrl.growFootprint, 
                                 badPixelMask, ctrl.usePixelWeights, ctrl.maxBadPixelFraction);
    } else {
        return ModelInputHandler(image, center, footprint, ctrl.growFootprint, 
                                 badPixelMask, ctrl.usePixelWeights, ctrl.maxBadPixelFraction);
    }
}

void FitProfileAlgorithm::fitShapeletTerms(
    FitProfileControl const & ctrl,
    FitPsfModel const & psfModel,
    ModelInputHandler const & inputs,
    FitProfileModel & model
) {
    typedef shapelet::MultiShapeletFunction MSF; 
    MSF msf = model.asMultiShapelet().convolve(psfModel.asMultiShapelet());
    msf.normalize();
    ndarray::Array<double,1,1> vector = ndarray::allocate(inputs.getSize());
    vector.deep() = 0.0;
    shapelet::ModelBuilder<double> builder(inputs.getX(), inputs.getY(), ctrl.useApproximateExp);
    for (MSF::ElementList::const_iterator i = msf.getElements().begin(); i != msf.getElements().end(); ++i) {
        builder.update(i->getEllipse().getCore());
        builder.addModelVector(i->getOrder(), i->getCoefficients(), vector);
    }
    if (!inputs.getWeights().isEmpty()) {
        vector.asEigen<Eigen::ArrayXpr>() *= inputs.getWeights().asEigen<Eigen::ArrayXpr>();
    }
    // the following is just linear least squares with one free parameter
    double variance = 1.0 / vector.asEigen().squaredNorm();
    model.flux = vector.asEigen().dot(inputs.getData().asEigen()) * variance;
    model.fluxErr = std::sqrt(variance);
    model.chisq = (model.flux * vector.asEigen() - inputs.getData().asEigen()).squaredNorm()
        / (inputs.getSize() - 4);
}

FitProfileModel FitProfileAlgorithm::apply(
    FitProfileControl const & ctrl,
    FitPsfModel const & psfModel,
    MultiGaussianObjective::EllipseCore const & inEllipse,
    ModelInputHandler const & inputs
) {
    HybridOptimizer opt = makeOptimizer(ctrl, psfModel, inEllipse, inputs);
    opt.run();
    Model model(
        ctrl, 
        boost::static_pointer_cast<MultiGaussianObjective const>(opt.getObjective())->getAmplitude(),
        opt.getParameters()
    );
    MultiGaussianObjective::EllipseCore ellipse = MultiGaussianObjective::readParameters(opt.getParameters());
    std::pair<bool,bool> constrained 
        = MultiGaussianObjective::constrainEllipse(ellipse, ctrl.minRadius, ctrl.minAxisRatio);
    model.flagMaxIter = opt.getState() & HybridOptimizer::FAILURE_MAXITER;
    model.flagTinyStep = (opt.getState() & HybridOptimizer::FAILURE_MINSTEP)
        || (opt.getState() & HybridOptimizer::FAILURE_MINTRUST);
    model.flagMinRadius = constrained.first;
    model.flagMinAxisRatio = constrained.second;
    model.flagLargeArea = !(model.ellipse.getArea() < inputs.getSize());
    model.fluxFlag = model.flagLargeArea;
    fitShapeletTerms(ctrl, psfModel, inputs, model);
    return model;
}

template <typename PixelT>
void FitProfileAlgorithm::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    source.set(_fluxKeys.flag, true);
    if (!exposure.hasPsf()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Cannot run FitProfileAlgorithm without a PSF."
        );
    }
    FitPsfModel psfModel(*_psfCtrl, source);
    if (psfModel.hasFailed() || !(psfModel.ellipse.getArea() > 0.0)) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "PSF shapelet fit failed; cannot fit galaxy model."
        );
    }
    afw::geom::ellipses::Quadrupole shape = source.getShape();
    if (source.getShapeFlag()) {
        shape = psfModel.ellipse;
    }
    ModelInputHandler inputs = adjustInputs(
        getControl(), psfModel, shape, *source.getFootprint(), exposure.getMaskedImage(), center
    );
    FitProfileModel model = apply(getControl(), psfModel, shape, inputs);

    assert(model.fluxFlag || lsst::utils::isfinite(model.ellipse.getArea()));

    source.set(_fluxKeys.meas, model.flux);
    source.set(_fluxKeys.err, model.fluxErr);
    source.set(_fluxKeys.flag, model.fluxFlag);
    source.set(_ellipseKey, model.ellipse);
    source.set(_chisqKey, model.chisq);
    source.set(_flagMaxIterKey, model.flagMaxIter);
    source.set(_flagTinyStepKey, model.flagTinyStep);
    source.set(_flagMinRadiusKey, model.flagMinRadius);
    source.set(_flagMinAxisRatioKey, model.flagMinAxisRatio);
    source.set(_flagLargeAreaKey, model.flagLargeArea);
}

template <typename PixelT>
void FitProfileAlgorithm::_applyForced(
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
            "Cannot run FitProfileAlgorithm without a PSF."
        );
    }
    if (source.getPsfFluxFlag()) {
        // The PSF likely has problems, which will cause us to chase our tail for a long time and use an
        // awful lot of memory in the process.  Better to just bail.
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "PSF flux failed; suspect galaxy modeling will be problematic."
        );
    }

    FitPsfModel psfModel(*_psfCtrl, source);
    if (psfModel.hasFailed() || !(psfModel.ellipse.getArea() > 0.0)) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "PSF shapelet fit failed; cannot fit galaxy model."
        );
    }

    FitProfileModel model(getControl(), reference);
    if (model.fluxFlag) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "Reference galaxy model fit failed; cannot run forced modeling."
        );
    }
    model.ellipse = model.ellipse.transform(refToMeas.getLinear());

    ModelInputHandler inputs = adjustInputs(
        getControl(), psfModel, model.ellipse, *source.getFootprint(),
        exposure.getMaskedImage(), center, true
    );

    fitShapeletTerms(getControl(), psfModel, inputs, model);

    source.set(_fluxKeys.meas, model.flux);
    source.set(_fluxKeys.err, model.fluxErr);
    source.set(_fluxKeys.flag, model.fluxFlag);
    source.set(_ellipseKey, model.ellipse);
    source.set(_chisqKey, model.chisq);
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(FitProfileAlgorithm);

}}}} // namespace lsst::meas::extensions::multiShapelet
