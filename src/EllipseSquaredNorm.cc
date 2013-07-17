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

#include "lsst/meas/extensions/multiShapelet/EllipseSquaredNorm.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

afw::geom::ellipses::BaseCore::Jacobian EllipseSquaredNorm::update(
    afw::geom::ellipses::BaseCore const & ellipse,
    bool computeJacobian
) {
    afw::geom::ellipses::Quadrupole q;
    afw::geom::ellipses::BaseCore::Jacobian dQdE = afw::geom::ellipses::BaseCore::Jacobian::Zero();
    if (computeJacobian) {
        dQdE = q.dAssign(ellipse);
    } else {
        q = ellipse;
    }
    // extract moments into variables just for readability
    double qxx = q.getIxx();
    double qxy = q.getIxy();
    double qyy = q.getIyy();
    // l = lower Cholesky decomposition of q.
    double l11 = std::sqrt(qxx);
    double l12 = qxy / l11;
    double l22 = std::sqrt(qyy - l12 * l12);
    // _r = inverse of l.
    _r11 = 1.0 / l11;
    _r22 = 1.0 / l22;
    _r12 = - _r11 * _r22 * l12;
    afw::geom::ellipses::BaseCore::Jacobian result;
    if (computeJacobian) {
        double s = l11 * l22; s *= s; // == det(Q)
        afw::geom::ellipses::BaseCore::Jacobian dRdQ;
        dRdQ << 
                          -0.5*_r11/qxx,             0.0,           0.0,
            -0.5*_r12*(1.0/qxx + qyy/s),  0.5*_r22*qxy/s,   -_r22*qyy/s,
                         0.5*_r12*qxy/s, -0.5*_r22*qxx/s,    _r22*qxy/s;
        result = dRdQ * dQdE;
    } else {
        result.setConstant(std::numeric_limits<double>::quiet_NaN());
    }
    return result;
}


}}}} // namespace lsst::meas::extensions::multiShapelet
