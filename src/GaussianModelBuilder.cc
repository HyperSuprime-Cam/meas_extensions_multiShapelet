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

#include "lsst/utils/PowFast.h"
#include "lsst/meas/extensions/multiShapelet/GaussianModelBuilder.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

namespace {

utils::PowFast const & powFast = utils::getPowFast<11>();

struct PowFastExpFunctor {
    inline float operator()(float x) const {
        return powFast.exp(x);
    }
};

} // anonymous

GaussianModelBuilder::GaussianModelBuilder(
    ndarray::Array<double const,1,1> const & x,
    ndarray::Array<double const,1,1> const & y,
    double flux, double radius, afw::geom::ellipses::Quadrupole const & psfEllipse,
    double psfAmplitude, bool useApproximateExp
) : _useApproximateExp(useApproximateExp), _flux(flux), _psfAmplitude(psfAmplitude), 
    _scaling(afw::geom::LinearTransform::makeScaling(radius)), _psfEllipse(psfEllipse), 
    _x(x), _y(y), _rx(x.size()), _ry(y.size())
{
    if (_x.size() != _y.size()) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("x coordinate array size (%d) does not match y coordinate array size (%d)")
             % _x.size() % _y.size()).str()
        );
    }
}

GaussianModelBuilder::GaussianModelBuilder(GaussianModelBuilder const & other) :
    _useApproximateExp(other._useApproximateExp), _flux(other._flux), _psfAmplitude(other._psfAmplitude),
    _scaling(other._scaling), _psfEllipse(other._psfEllipse),
    _esn(other._esn), _esnJacobian(other._esnJacobian), _dNorm(other._dNorm),
    _x(other._x), _y(other._y), _rx(other._rx), _ry(other._ry)
{}

GaussianModelBuilder & GaussianModelBuilder::operator=(GaussianModelBuilder const & other) {
    if (&other != this) {
        _useApproximateExp = other._useApproximateExp;
        _flux = other._flux;
        _psfAmplitude = other._psfAmplitude;
        _scaling = other._scaling;
        _psfEllipse = other._psfEllipse;
        _x.reset(other._x.shallow());
        _y.reset(other._y.shallow());
        _rx = other._rx;
        _ry = other._ry;
        _esn = other._esn;
        _esnJacobian = other._esnJacobian;
        _dNorm = other._dNorm;
        if (!other._model.isEmpty()) _model = ndarray::copy(other._model);
    }
    return *this;
}

void GaussianModelBuilder::update(afw::geom::ellipses::BaseCore const & core) {
    Eigen::Matrix3d scaleJac = core.transform(_scaling).d();
    PTR(afw::geom::ellipses::BaseCore) ellipse = core.transform(_scaling).copy();
    Eigen::Matrix3d convJac = ellipse->convolve(_psfEllipse).d();
    ellipse->convolve(_psfEllipse).inPlace();
    afw::geom::ellipses::Quadrupole q;
    Eigen::Matrix3d quadJacobian = q.dAssign(*ellipse) * convJac * scaleJac;
    _esnJacobian = _esn.update(q) * quadJacobian; 
    if (_model.isEmpty()) {
        _model = ndarray::allocate(_x.size());
    }
    ndarray::EigenView<double,1,1> z(_model);
    _esn(_x, _y, _rx, _ry, z);
    double det = q.getDeterminant();
    Eigen::RowVector3d dNorm_dq;
    dNorm_dq[0] = -0.5 * q.getIyy() / det;
    dNorm_dq[1] = -0.5 * q.getIxx() / det;
    dNorm_dq[2] = q.getIxy() / det;
    _dNorm = dNorm_dq * quadJacobian;
    if (_useApproximateExp) {
        z.array() = (-0.5 * z.array()).unaryExpr(PowFastExpFunctor());
    } else {
        z.array() = (-0.5 * z.array()).exp();
    }
    z.array() *= (_flux * _psfAmplitude) / (std::sqrt(det) * afw::geom::PI * 2.0);
}

void GaussianModelBuilder::computeDerivative(
    ndarray::Array<double,2,-1> const & output,
    bool add
) {
    if (_model.isEmpty()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "computeDerivative called before computeModel"
        );
    }
    if (output.getSize<0>() != _x.size()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            (boost::format("Incorrect number of rows for array: got %d, expected %d")
             % output.getSize<0>() % _x.size()).str()
        );
    }
    ndarray::EigenView<double,2,-1> out(output);
    if (!add) out.setZero();
    Eigen::MatrixXd dz_de = Eigen::MatrixXd::Zero(_x.size(), _esnJacobian.cols());
    _esn.dEllipse(_x, _y, _rx, _ry, _esnJacobian, dz_de);
    for (int n = 0; n < _esnJacobian.cols(); ++n) {
        dz_de.col(n).array() *= -0.5 * _model.asEigen<Eigen::ArrayXpr>();
    }
    out += dz_de;
    out += _model.asEigen() * _dNorm;
}

void GaussianModelBuilder::setOutput(ndarray::Array<double,1,1> const & array) {
    if (array.getSize<0>() != _x.size()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            (boost::format("Incorrect size for array: got %d, expected %d")
             % array.getSize<0>() % _x.size()).str()
        );
    }
}

}}}} // namespace lsst::meas::extensions::multiShapelet
