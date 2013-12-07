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

#include "lsst/meas/extensions/multiShapelet/FitPsf.h"
#include "lsst/meas/extensions/multiShapelet/MultiGaussianObjective.h"
#include "lsst/shapelet/ModelBuilder.h"
#include "lsst/afw/math/LeastSquares.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

namespace {

int computeOrder(int size2d) {
    // Could just use f.p. square roots here, but we only care about small numbers
    // and I don't want to think about how to be robust against round-off error issues.
    for (int o=0, s=1; s <= size2d; s += (++o + 1)) {
        if (s == size2d) return o;
    }
    throw LSST_EXCEPT(
        pex::exceptions::LogicErrorException,
        (boost::format("%d is not a valid shapelet vector size") % size2d).str()
    );
}

} // anonymous

MultiGaussian FitPsfControl::getMultiGaussian() const {
    MultiGaussian multiGaussian;
    multiGaussian.add(GaussianComponent(1.0, 1.0));
    multiGaussian.add(GaussianComponent(peakRatio * radiusRatio * radiusRatio, radiusRatio));
    return multiGaussian;
}

FitPsfModel::FitPsfModel(
    FitPsfControl const & ctrl,
    double amplitude,
    ndarray::Array<double const,1,1> const & parameters
) :
    ellipse(),
    radiusRatio(ctrl.radiusRatio), chisq(std::numeric_limits<double>::quiet_NaN()),
    failedMaxIter(false), failedTinyStep(false), failedMinRadius(false), failedMinAxisRatio(false)
{
    MultiGaussian multiGaussian = ctrl.getMultiGaussian();
    ellipse = MultiGaussianObjective::EllipseCore(parameters[0], parameters[1], parameters[2]);
    shapelet::ShapeletFunction innerShapelet = multiGaussian[0].makeShapelet(
        afw::geom::ellipses::Ellipse(ellipse), ctrl.innerOrder
    );
    inner = innerShapelet.getCoefficients();
    inner.asEigen() *= amplitude;
    shapelet::ShapeletFunction outerShapelet = multiGaussian[1].makeShapelet(
        afw::geom::ellipses::Ellipse(ellipse), ctrl.outerOrder
    );
    outer = outerShapelet.getCoefficients();
    outer.asEigen() *= amplitude;
}

FitPsfModel::FitPsfModel(FitPsfControl const & ctrl, afw::table::BaseRecord const & source) :
    inner(ndarray::allocate(shapelet::computeSize(ctrl.innerOrder))),
    outer(ndarray::allocate(shapelet::computeSize(ctrl.outerOrder))),
    ellipse(),
    radiusRatio(ctrl.radiusRatio), chisq(std::numeric_limits<double>::quiet_NaN()),
    failedMaxIter(false), failedTinyStep(false), failedMinRadius(false), failedMinAxisRatio(false)
{
    afw::table::SubSchema s = source.getSchema()[ctrl.name];
    afw::table::Key< afw::table::Array<float> > innerKey = s.find< afw::table::Array<float> >("inner").key;
    if (innerKey.getSize() <= inner.getSize<0>()) {
        inner.deep() = 0.0;
        inner[ndarray::view(0, innerKey.getSize())] = source.get(innerKey);
    } else {
        inner.deep() = source.get(innerKey)[ndarray::view(0, inner.getSize<0>())];
    }
    afw::table::Key< afw::table::Array<float> > outerKey = s.find< afw::table::Array<float> >("outer").key; 
    if (outerKey.getSize() <= outer.getSize<0>()) {
        outer.deep() = 0.0;
        outer[ndarray::view(0, outerKey.getSize())] = source.get(outerKey);
    } else {
        outer.deep() = source.get(outerKey)[ndarray::view(0, outer.getSize<0>())];
    }
    chisq = source.get(s.find< float >("chisq").key);
    ellipse = source.get(s.find< afw::table::Moments<float> >("ellipse").key);
    failedMaxIter = source.get(s.find<afw::table::Flag>("flags.maxiter").key);
    failedTinyStep = source.get(s.find<afw::table::Flag>("flags.tinystep").key);
    failedMinRadius = source.get(s.find<afw::table::Flag>("flags.constraint.r").key);
    failedMinRadius = source.get(s.find<afw::table::Flag>("flags.constraint.q").key);
}

FitPsfModel::FitPsfModel(FitPsfModel const & other) :
    inner(ndarray::copy(other.inner)),
    outer(ndarray::copy(other.outer)),
    ellipse(other.ellipse),
    radiusRatio(other.radiusRatio),
    chisq(other.chisq),
    failedMaxIter(other.failedMaxIter),
    failedTinyStep(other.failedTinyStep),
    failedMinRadius(other.failedMinRadius),
    failedMinAxisRatio(other.failedMinAxisRatio)    
{}

FitPsfModel & FitPsfModel::operator=(FitPsfModel const & other) {
    if (&other != this) {
        inner = ndarray::copy(other.inner);
        outer = ndarray::copy(other.outer);
        ellipse = other.ellipse;
        radiusRatio = other.radiusRatio;
        chisq = other.chisq;
        failedMaxIter = other.failedMaxIter;
        failedTinyStep = other.failedTinyStep;
        failedMinRadius = other.failedMinRadius;
        failedMinAxisRatio = other.failedMinAxisRatio;
    }
    return *this;
}

MultiGaussian FitPsfModel::getMultiGaussian() const {
    MultiGaussian result;
    result.add(GaussianComponent(1.0, 1.0)).readShapeletAmplitude(inner[0], ellipse);
    result.add(GaussianComponent(1.0, radiusRatio)).readShapeletAmplitude(outer[0], ellipse);
    return result;
}

shapelet::MultiShapeletFunction FitPsfModel::asMultiShapelet(
    afw::geom::Point2D const & center
) const {
    shapelet::MultiShapeletFunction::ElementList elements;
    afw::geom::ellipses::Ellipse fullEllipse(ellipse, center);
    elements.push_back(
        shapelet::ShapeletFunction(
            computeOrder(inner.getSize<0>()),
            shapelet::HERMITE,
            fullEllipse,
            inner
        )
    );
    fullEllipse.scale(radiusRatio);
    elements.push_back(
        shapelet::ShapeletFunction(
            computeOrder(outer.getSize<0>()),
            shapelet::HERMITE,
            fullEllipse,
            outer
        )
    );
    return shapelet::MultiShapeletFunction(elements);
}

FitPsfAlgorithm::FitPsfAlgorithm(FitPsfControl const & ctrl, afw::table::Schema & schema) :
    algorithms::Algorithm(ctrl),
    _innerKey(
        schema.addField< afw::table::Array<float> >(
            ctrl.name + ".inner",
            "Gauss-Hermite coefficients of the inner expansion (see lsst.shapelet package)",
            shapelet::computeSize(ctrl.innerOrder)
        )),
    _outerKey(
        schema.addField< afw::table::Array<float> >(
            ctrl.name + ".outer",
            "Gauss-Hermite coefficients of the outer expansion (see lsst.shapelet package)",
            shapelet::computeSize(ctrl.outerOrder)
        )),
    _ellipseKey(
        schema.addField< afw::table::Moments<float> >(
            ctrl.name + ".ellipse",
            "Ellipse corresponding to the inner expansion"
        )),
    _chisqKey(
        schema.addField<float>(ctrl.name + ".chisq", "Reduced chi^2 of the final shapelet fit")
    ),
    _integralKey(
        schema.addField<float>(
            ctrl.name + ".integral",
            "Integral of the shapelet PSF model to infinite radius"
        )
    ),
    _flagKey(
        schema.addField<afw::table::Flag>(
            ctrl.name + ".flags",
            "set if the multi-shapelet PSF fit was unsuccessful in any way"
        )),
    _flagMaxIterKey(
        schema.addField<afw::table::Flag>(
            ctrl.name + ".flags.maxiter",
            "set if the optimizer ran into the maximum number of iterations limit"
        )),
    _flagTinyStepKey(
        schema.addField<afw::table::Flag>(
            ctrl.name + ".flags.tinystep",
            "set if the optimizer step or trust region got so small no progress could be made"
        )),
    _flagMinRadiusKey(
        schema.addField<afw::table::Flag>(
            ctrl.name + ".flags.constraint.r",
            "set if the best-fit radius was the minimum allowed by the constraint"
        )),
    _flagMinAxisRatioKey(
        schema.addField<afw::table::Flag>(
            ctrl.name + ".flags.constraint.q",
            "set if the best-fit axis ratio (b/a) was the minimum allowed by the constraint"
        ))
{}

PTR(MultiGaussianObjective) FitPsfAlgorithm::makeObjective(
    FitPsfControl const & ctrl,
    ModelInputHandler const & inputs
) {
    return boost::make_shared<MultiGaussianObjective>(
        inputs, ctrl.getMultiGaussian(), ctrl.minRadius, ctrl.minAxisRatio, ctrl.useApproximateExp
    );
}

HybridOptimizer FitPsfAlgorithm::makeOptimizer(
    FitPsfControl const & ctrl,
    ModelInputHandler const & inputs
) {
    MultiGaussianObjective::EllipseCore ellipse(0.0, 0.0, std::log(ctrl.initialRadius));
    MultiGaussianObjective::constrainEllipse(ellipse, ctrl.minRadius, ctrl.minAxisRatio);
    PTR(Objective) obj = makeObjective(ctrl, inputs);
    HybridOptimizerControl optCtrl; // TODO: nest this in FitPsfControl
    optCtrl.tau = 1E-6;
    optCtrl.useCholesky = true;
    optCtrl.gTol = 1E-6;
    ndarray::Array<double,1,1> initial = ndarray::allocate(obj->getParameterSize());
    ellipse.writeParameters(initial.getData());
    return HybridOptimizer(obj, initial, optCtrl);
}

void FitPsfAlgorithm::fitShapeletTerms(
    FitPsfControl const & ctrl,
    ModelInputHandler const & inputs,
    FitPsfModel & model
) {
    int innerCoeffs = shapelet::computeSize(ctrl.innerOrder);
    int outerCoeffs = shapelet::computeSize(ctrl.outerOrder);
    ndarray::Array<double,2,-2> matrix = ndarray::allocate(inputs.getSize(), innerCoeffs + outerCoeffs);
    matrix.asEigen().setZero();
    shapelet::ModelBuilder<double> builder(inputs.getX(), inputs.getY(), ctrl.useApproximateExp);
    builder.update(model.ellipse);
    builder.addModelMatrix(ctrl.innerOrder, matrix[ndarray::view()(0, innerCoeffs)]);
    model.ellipse.scale(ctrl.radiusRatio);
    builder.update(model.ellipse);
    model.ellipse.scale(1.0 / ctrl.radiusRatio);
    builder.addModelMatrix(ctrl.outerOrder, matrix[ndarray::view()(innerCoeffs, innerCoeffs + outerCoeffs)]);
    if (!inputs.getWeights().isEmpty()) {
        matrix.asEigen<Eigen::ArrayXpr>() 
            *= (inputs.getWeights().asEigen() * Eigen::RowVectorXd::Ones(matrix.getSize<1>())).array();
    }
    afw::math::LeastSquares lstsq = afw::math::LeastSquares::fromDesignMatrix(matrix, inputs.getData());
    model.inner.deep() = lstsq.getSolution()[ndarray::view(0, innerCoeffs)];
    model.outer.deep() = lstsq.getSolution()[ndarray::view(innerCoeffs, innerCoeffs + outerCoeffs)];
    // The degrees of freedom corresponds to the final shapelet fit with ellipse held fixed.
    model.chisq =
        (matrix.asEigen() * lstsq.getSolution().asEigen() - inputs.getData().asEigen()).squaredNorm()
        / (matrix.getSize<0>() - matrix.getSize<1>());
}

FitPsfModel FitPsfAlgorithm::apply(
    FitPsfControl const & ctrl,
    ModelInputHandler const & inputs
) {
    HybridOptimizer opt = makeOptimizer(ctrl, inputs);
    opt.run();
    Model model(
        ctrl, 
        boost::static_pointer_cast<MultiGaussianObjective const>(opt.getObjective())->getAmplitude(),
        opt.getParameters()
    );
    MultiGaussianObjective::EllipseCore ellipse = MultiGaussianObjective::readParameters(opt.getParameters());
    std::pair<bool,bool> constrained 
        = MultiGaussianObjective::constrainEllipse(ellipse, ctrl.minRadius, ctrl.minAxisRatio);
    model.failedMaxIter = opt.getState() & HybridOptimizer::FAILURE_MAXITER;
    model.failedTinyStep = (opt.getState() & HybridOptimizer::FAILURE_MINSTEP)
        || (opt.getState() & HybridOptimizer::FAILURE_MINTRUST);
    model.failedMinRadius = constrained.first;
    model.failedMinAxisRatio = constrained.second;
    fitShapeletTerms(ctrl, inputs, model);
    return model;
}

FitPsfModel FitPsfAlgorithm::apply(
    FitPsfControl const & ctrl,
    afw::detection::Psf const & psf,
    afw::geom::Point2D const & center
) {
    PTR(afw::image::Image<afw::math::Kernel::Pixel>) image = psf.computeImage(center);
    double s = image->getArray().asEigen().sum();
    image->getArray().asEigen() /= s;
    ModelInputHandler inputs(*image, center, image->getBBox(afw::image::PARENT));
    return apply(ctrl, inputs);
}

FitPsfModel FitPsfAlgorithm::fit(
    afw::table::BaseRecord & record,
    afw::detection::Psf const & psf,
    afw::geom::Point2D const & center
) const {
    record.set(_flagKey, true);
    FitPsfModel model = apply(getControl(), psf, center);
    record[_innerKey] = model.inner;
    record[_outerKey] = model.outer;
    record.set(_ellipseKey, model.ellipse);
    record.set(_chisqKey, model.chisq);
    record.set(_integralKey, model.asMultiShapelet().evaluate().integrate());
    record.set(_flagMaxIterKey, model.failedMaxIter);
    record.set(_flagTinyStepKey, model.failedTinyStep);
    record.set(_flagMinRadiusKey, model.failedMinRadius);
    record.set(_flagMinAxisRatioKey, model.failedMinAxisRatio);
    record.set(_flagKey, model.failedMaxIter || model.failedTinyStep
               || model.failedMinAxisRatio || model.failedMinRadius);
    return model;
}

template <typename PixelT>
void FitPsfAlgorithm::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    source.set(_flagKey, true);
    if (!exposure.hasPsf()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Cannot run FitPsfAlgorithm without a PSF."
        );
    }
    fit(source, *exposure.getPsf(), center);
}

PTR(algorithms::AlgorithmControl) FitPsfControl::_clone() const {
    return boost::make_shared<FitPsfControl>(*this);
}

PTR(algorithms::Algorithm) FitPsfControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::make_shared<FitPsfAlgorithm>(*this, boost::ref(schema));
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(FitPsfAlgorithm);


}}}} // namespace lsst::meas::extensions::multiShapelet
