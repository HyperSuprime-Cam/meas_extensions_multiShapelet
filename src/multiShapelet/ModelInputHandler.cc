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

#include "ndarray/eigen.h"
#include "lsst/meas/extensions/multiShapelet/ModelInputHandler.h"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

namespace {

void initCoords(
    ndarray::Array<double,1,1> & x, ndarray::Array<double,1,1> & y,
    afw::detection::Footprint const & region, afw::geom::Point2D const & center
) {
    x = ndarray::allocate(region.getArea());
    y = ndarray::allocate(region.getArea());
    int n = 0;
    for (
        afw::detection::Footprint::SpanList::const_iterator spanIter = region.getSpans().begin();
        spanIter != region.getSpans().end();
        ++spanIter
    ) {
        double dy = (**spanIter).getY() - center.getY();
        for (int ix = (**spanIter).getX0(); ix <= (**spanIter).getX1(); ++ix) {
            x[n] = ix - center.getX();
            y[n] = dy;
            ++n;
        }
    }
}

} // anonymous

template <typename PixelT>
ModelInputHandler::ModelInputHandler(
    afw::image::Image<PixelT> const & image, afw::geom::Point2D const & center, 
    afw::geom::Box2I const & region
) {
    _footprint = boost::make_shared<afw::detection::Footprint>(region);
    _footprint->clipTo(image.getBBox(afw::image::PARENT));
    _data = ndarray::allocate(_footprint->getArea());
    afw::detection::flattenArray(*_footprint, image.getArray(), _data, image.getXY0());
    initCoords(_x, _y, *_footprint, center);
}

template <typename PixelT>
ModelInputHandler::ModelInputHandler(
    afw::image::Image<PixelT> const & image, afw::geom::Point2D const & center, 
    afw::detection::Footprint const & region, int growFootprint
) {
    if (growFootprint) {
        _footprint = afw::detection::growFootprint(region, growFootprint);
    } else {
        _footprint = boost::make_shared<afw::detection::Footprint>(region);
    }
    _footprint->clipTo(image.getBBox(afw::image::PARENT));
    _data = ndarray::allocate(_footprint->getArea());
    afw::detection::flattenArray(*_footprint, image.getArray(), _data, image.getXY0());
    initCoords(_x, _y, *_footprint, center);
}

template <typename PixelT>
ModelInputHandler::ModelInputHandler(
    afw::image::MaskedImage<PixelT> const & image, afw::geom::Point2D const & center, 
    afw::geom::Box2I const & region, afw::image::MaskPixel badPixelMask, bool usePixelWeights
) {
    _footprint = boost::make_shared<afw::detection::Footprint>(region);
    _footprint->intersectMask(*image.getMask(), badPixelMask);
    _data = ndarray::allocate(_footprint->getArea());
    _weights = ndarray::allocate(_footprint->getArea());
    afw::detection::flattenArray(*_footprint, image.getImage()->getArray(), _data, image.getXY0());
    afw::detection::flattenArray(*_footprint, image.getVariance()->getArray(), _weights, image.getXY0());
    if (!usePixelWeights) {
        _weights.asEigen().setConstant(_weights.asEigen().mean());
    }
    _weights.asEigen<Eigen::ArrayXpr>() = _weights.asEigen<Eigen::ArrayXpr>().sqrt().inverse();
    _data.asEigen<Eigen::ArrayXpr>() *= _weights.asEigen<Eigen::ArrayXpr>();
    initCoords(_x, _y, *_footprint, center);
}

template <typename PixelT>
ModelInputHandler::ModelInputHandler(
    afw::image::MaskedImage<PixelT> const & image, afw::geom::Point2D const & center, 
    afw::detection::Footprint const & region, int growFootprint,
    afw::image::MaskPixel badPixelMask, bool usePixelWeights
) {
    if (growFootprint) {
        _footprint = afw::detection::growFootprint(region, growFootprint);
    } else {
        _footprint = boost::make_shared<afw::detection::Footprint>(region);
    }
    _footprint->intersectMask(*image.getMask(), badPixelMask);
    _data = ndarray::allocate(_footprint->getArea());
    _weights = ndarray::allocate(_footprint->getArea());
    afw::detection::flattenArray(*_footprint, image.getImage()->getArray(), _data, image.getXY0());
    afw::detection::flattenArray(*_footprint, image.getVariance()->getArray(), _weights, image.getXY0());
    if (!usePixelWeights) {
        _weights.asEigen().setConstant(_weights.asEigen().mean());
    }
    _weights.asEigen<Eigen::ArrayXpr>() = _weights.asEigen<Eigen::ArrayXpr>().sqrt().inverse();
    _data.asEigen<Eigen::ArrayXpr>() *= _weights.asEigen<Eigen::ArrayXpr>();
    initCoords(_x, _y, *_footprint, center);
}

#define INSTANTIATE(T)                          \
    template ModelInputHandler::ModelInputHandler(                      \
        afw::image::Image<T> const & image, afw::geom::Point2D const & center, \
        afw::geom::Box2I const & region);                               \
    template ModelInputHandler::ModelInputHandler(                      \
        afw::image::Image<T> const & image, afw::geom::Point2D const & center, \
        afw::detection::Footprint const & region, int growFootprint);   \
    template ModelInputHandler::ModelInputHandler(                      \
        afw::image::MaskedImage<T> const & image, afw::geom::Point2D const & center, \
        afw::geom::Box2I const & region, afw::image::MaskPixel badPixelMask, bool usePixelWeights); \
    template ModelInputHandler::ModelInputHandler(                      \
        afw::image::MaskedImage<T> const & image, afw::geom::Point2D const & center, \
        afw::detection::Footprint const & region, int growFootprint,  \
        afw::image::MaskPixel badPixelMask, bool usePixelWeights)

INSTANTIATE(float);
INSTANTIATE(double);

}}}} // namespace lsst::meas::extensions::multiShapelet
