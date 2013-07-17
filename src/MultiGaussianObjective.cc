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
#include "lsst/meas/extensions/multiShapelet/MultiGaussianObjective.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

MultiGaussianObjective::MultiGaussianObjective(
    ModelInputHandler const & inputs,
    MultiGaussian const & multiGaussian,
    double minRadius, double minAxisRatio,
    bool useApproximateExp
) : Objective(inputs.getSize(), 3), _minRadius(minRadius), _minAxisRatio(minAxisRatio), 
    _amplitude(1.0), _modelSquaredNorm(1.0),
    _ellipse(), _inputs(inputs), _model(ndarray::allocate(inputs.getSize()))
{
    if (_minRadius <= 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            "Minimum radius must be = 0"
        );
    }
    if (_minAxisRatio < 0.0 || _minAxisRatio > 1.0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            "Minimum axis ratio must be between 0 and 1"
        );
    }
    _builders.reserve(multiGaussian.size());
    for (MultiGaussian::const_iterator i = multiGaussian.begin(); i != multiGaussian.end(); ++i) {
        _builders.push_back(
            GaussianModelBuilder(
                _inputs.getX(), _inputs.getY(), i->flux, i->radius,
                afw::geom::ellipses::Quadrupole(0.0, 0.0, 0.0), 1.0,
                useApproximateExp
            )
        );
    }
}

MultiGaussianObjective::MultiGaussianObjective(
    ModelInputHandler const & inputs,
    MultiGaussian const & multiGaussian,
    MultiGaussian const & psfMultiGaussian,
    afw::geom::ellipses::Quadrupole const & psfEllipse,
    double minRadius, double minAxisRatio,
    bool useApproximateExp
) : Objective(inputs.getSize(), 3), _minRadius(minRadius), _minAxisRatio(minAxisRatio),
    _amplitude(1.0), _modelSquaredNorm(1.0),
    _ellipse(), _inputs(inputs), _model(ndarray::allocate(inputs.getSize()))
{
    if (_minRadius <= 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            "Minimum radius must be > 0"
        );
    }
    if (_minAxisRatio < 0.0 || _minAxisRatio > 1.0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            "Minimum axis ratio must be between 0 and 1"
        );
    }
    _builders.reserve(multiGaussian.size() * psfMultiGaussian.size());
    for (MultiGaussian::const_iterator j = psfMultiGaussian.begin(); j != psfMultiGaussian.end(); ++j) {
        afw::geom::ellipses::Quadrupole psfComponentEllipse(psfEllipse);
        psfComponentEllipse.scale(j->radius);
        for (MultiGaussian::const_iterator i = multiGaussian.begin(); i != multiGaussian.end(); ++i) {
            _builders.push_back(
                GaussianModelBuilder(
                    _inputs.getX(), _inputs.getY(), i->flux, i->radius,
                    psfComponentEllipse, j->flux,
                    useApproximateExp
                )
            );
        }
    }
}

Objective::StepResult MultiGaussianObjective::tryStep(
    ndarray::Array<double const,1,1> const & oldParameters, 
    ndarray::Array<double,1,1> const & newParameters
) {
    StepResult result(VALID);
    for (int n = 0; n < oldParameters.getSize<0>(); ++n) {
        if (!lsst::utils::isfinite(newParameters[n]) && lsst::utils::isfinite(oldParameters[n])) {
            newParameters[n] = oldParameters[n];
        }
    }
    _ellipse.readParameters(newParameters.getData());
    std::pair<bool,bool> constrained = constrainEllipse(_ellipse, _minRadius, _minAxisRatio);
    if (constrained.first || constrained.second) {
        result = MODIFIED;
        _ellipse.writeParameters(newParameters.getData());
    }
    return result;
}

void MultiGaussianObjective::computeFunction(
    ndarray::Array<double const,1,1> const & parameters, 
    ndarray::Array<double,1,1> const & function
) {
    _ellipse.readParameters(parameters.getData());
    ndarray::EigenView<double,1,1> model(_model);
    model.setZero();
    for (std::size_t n = 0; n < _builders.size(); ++n) {
        _builders[n].update(_ellipse);
        model += _builders[n].getModel().asEigen();
    }
    if (!_inputs.getWeights().isEmpty()) {
        model.array() *= _inputs.getWeights().asEigen<Eigen::ArrayXpr>();
    }
    _modelSquaredNorm = model.squaredNorm();
    _amplitude = model.dot(_inputs.getData().asEigen()) / _modelSquaredNorm;
    function.asEigen() = _amplitude * model - _inputs.getData().asEigen();
}

void MultiGaussianObjective::computeDerivative(
    ndarray::Array<double const,1,1> const & parameters, 
    ndarray::Array<double const,1,1> const & function,
    ndarray::Array<double,2,-2> const & derivative
) {
    derivative.asEigen().setZero();
    for (std::size_t n = 0; n < _builders.size(); ++n) {
        _builders[n].computeDerivative(derivative, true);
    }
    if (!_inputs.getWeights().isEmpty()) {
        derivative.asEigen<Eigen::ArrayXpr>() 
            *= (_inputs.getWeights().asEigen() * Eigen::RowVectorXd::Ones(parameters.getSize<0>())).array();
    }
    // Right now, 'derivative' is the partial derivative w.r.t. the objective parameters
    // with flux held fixed at 1.  However, the parameters also affect the flux, so we need
    // to compute the partial derivative of that.
    Eigen::VectorXd tmp = _inputs.getData().asEigen() - 2.0 * _amplitude * _model.asEigen();
    Eigen::VectorXd dAmplitude = (derivative.asEigen().adjoint() * tmp) / _modelSquaredNorm;
    // Now we update 'derivative' so it becomes the complete derivative rather than the partial.
    derivative.asEigen() *= _amplitude;
    derivative.asEigen() += _model.asEigen() * dAmplitude.transpose();
}


MultiGaussianObjective::EllipseCore
MultiGaussianObjective::readParameters(ndarray::Array<double const,1,1> const & parameters) {
    EllipseCore r;
    r.readParameters(parameters.getData());
    return r;
}

void MultiGaussianObjective::writeParameters(
    EllipseCore const & ellipse,
    ndarray::Array<double,1,1> const & parameters
) {
    ellipse.writeParameters(parameters.getData());
}

std::pair<bool,bool> 
MultiGaussianObjective::constrainEllipse(EllipseCore & ellipse, double minRadius, double minAxisRatio) {
    if (minRadius <= 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            "Minimum radius must be > 0"
        );
    }
    if (minAxisRatio < 0.0 || minAxisRatio > 1.0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            "Minimum axis ratio must be between 0 and 1"
        );
    }
    std::pair<bool,bool> result(false, false);
    // First we make sure we don't have any NaNs or Infs; using ellipse setters with those
    // present may cause non-finite values to spread to other parameters.
    double e1 = ellipse.getEllipticity().getE1();
    double e2 = ellipse.getEllipticity().getE2();
    double r = ellipse.getRadius();
    if (!lsst::utils::isfinite(e1)) {
        result.second = true;
        e1 = 0.0;
    }
    if (!lsst::utils::isfinite(e2)) {
        result.second = true;
        e2 = 0.0;
    }
    if (!lsst::utils::isfinite(r)) {
        result.first = true;
        r = std::log(minRadius);
    }
    if (result.first || result.second) {
        ellipse = EllipseCore(e1, e2, r);
    }
    // Now we can move on to comparing with the constraints.  Note that we transform
    // the constraints into the logarithmic definition used by the parameters rather
    // than the reverse.
    // We use epsilon below to make sure we return true when the value is already at the constraint,
    // even if it may have been through some ellipse reparameterizations since then.
    if (minAxisRatio > std::numeric_limits<double>::epsilon()) {
        double eMax = -std::log(minAxisRatio);
        if (ellipse.getEllipticity().getE() >= eMax - std::numeric_limits<double>::epsilon()) {
            ellipse.getEllipticity().setE(eMax);
            result.second = true;
        }
    }
    if (minRadius > std::numeric_limits<double>::epsilon()) {
        double rMin = std::log(minRadius);
        if (ellipse.getRadius() <= rMin + std::numeric_limits<double>::epsilon()) {
            ellipse.setRadius(rMin);
            result.first = true;
        }
    }
    return result;
}

}}}} // namespace lsst::meas::extensions::multiShapelet
