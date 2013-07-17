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
#include "lsst/afw/detection/FootprintSet.h"
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

PTR(afw::detection::Footprint) mergeFootprintWithEllipses(
    afw::detection::Footprint const & footprint,
    std::vector<afw::geom::ellipses::Ellipse> const & ellipses,
    afw::geom::Box2I const & imageBox
) {
    // TODO: we do a lot of turning footprints into masks here and elsewhere; should really only
    // have to do that once
    afw::geom::Box2I bbox(footprint.getBBox());
    std::vector<PTR(afw::detection::Footprint)> ellipseFootprints(ellipses.size());
    for (std::size_t n = 0; n < ellipses.size(); ++n) {
        ellipseFootprints[n] = boost::make_shared<afw::detection::Footprint>(ellipses[n]);
        bbox.include(ellipseFootprints[n]->getBBox());
    }

    if (!(bbox.getArea() > 0) || !(imageBox.contains(bbox))) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "Invalid bounding box in model fit"
        );
    }
    afw::image::Mask<> mask(bbox);
    assert(mask.getBBox(afw::image::PARENT).contains(footprint.getBBox()));
    afw::detection::setMaskFromFootprint(&mask, footprint, afw::image::MaskPixel(0x1));
    for (std::size_t n = 0; n < ellipses.size(); ++n) {
        assert(mask.getBBox(afw::image::PARENT).contains(ellipseFootprints[n]->getBBox()));
        afw::detection::setMaskFromFootprint(&mask, *ellipseFootprints[n], afw::image::MaskPixel(0x1));
    }
    afw::detection::FootprintSet fpSet(
        mask, afw::detection::Threshold(0x1, afw::detection::Threshold::BITMASK), 1
    );
    if (fpSet.getFootprints()->size() != 1u) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            "Ellipse-based footprints do not all overlap detection footprint."
        );
    }
    return fpSet.getFootprints()->front();
}

} // anonymous

template <typename PixelT>
ModelInputHandler::ModelInputHandler(
    afw::image::Image<PixelT> const & image, afw::geom::Point2D const & center, 
    afw::geom::Box2I const & region
) {
    _footprint = boost::make_shared<afw::detection::Footprint>(region);
    _footprint->clipTo(image.getBBox(afw::image::PARENT));
    if (_footprint->getArea() <= 0) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "Model fit contains no usable pixels."
        );
    }
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
    if (_footprint->getArea() <= 0) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "Model fit contains no usable pixels."
        );
    }
    _data = ndarray::allocate(_footprint->getArea());
    afw::detection::flattenArray(*_footprint, image.getArray(), _data, image.getXY0());
    initCoords(_x, _y, *_footprint, center);
}

template <typename PixelT>
ModelInputHandler::ModelInputHandler(
    afw::image::Image<PixelT> const & image, afw::geom::Point2D const & center,
    std::vector<afw::geom::ellipses::Ellipse> const & ellipses, 
    afw::detection::Footprint const & region, int growFootprint
) {
    if (growFootprint) {
        _footprint = afw::detection::growFootprint(region, growFootprint);
    } else {
        _footprint = boost::make_shared<afw::detection::Footprint>(region);
    }
    _footprint = mergeFootprintWithEllipses(*_footprint, ellipses, image.getBBox(afw::image::PARENT));
    _footprint->clipTo(image.getBBox(afw::image::PARENT));
    if (_footprint->getArea() <= 0) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "Model fit contains no usable pixels."
        );
    }
    _data = ndarray::allocate(_footprint->getArea());
    afw::detection::flattenArray(*_footprint, image.getArray(), _data, image.getXY0());
    initCoords(_x, _y, *_footprint, center);
}

template <typename PixelT>
ModelInputHandler::ModelInputHandler(
    afw::image::MaskedImage<PixelT> const & image, afw::geom::Point2D const & center, 
    afw::geom::Box2I const & region, afw::image::MaskPixel badPixelMask, bool usePixelWeights,
    double maxBadPixelFraction
) {
    _footprint = boost::make_shared<afw::detection::Footprint>(region);
    double originalArea = _footprint->getArea();
    _footprint->intersectMask(*image.getMask(), badPixelMask);
    if ((1.0 - _footprint->getArea() / originalArea) > maxBadPixelFraction) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "Too many masked pixels."
        );
    } 
    if (_footprint->getArea() <= 0) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "Model fit contains no usable pixels."
        );
    }
    _data = ndarray::allocate(_footprint->getArea());
    _weights = ndarray::allocate(_footprint->getArea());
    afw::detection::flattenArray(*_footprint, image.getImage()->getArray(), _data, image.getXY0());
    afw::detection::flattenArray(*_footprint, image.getVariance()->getArray(), _weights, image.getXY0());
    if (!usePixelWeights) {
        _weights.asEigen().setConstant(_weights.asEigen().mean());
    }
    // "operator=" needed here to workaround clang bug in resolving inherited assignment operators
    _weights.asEigen<Eigen::ArrayXpr>().operator=(_weights.asEigen<Eigen::ArrayXpr>().sqrt().inverse());
    _data.asEigen<Eigen::ArrayXpr>() *= _weights.asEigen<Eigen::ArrayXpr>();
    initCoords(_x, _y, *_footprint, center);
}

template <typename PixelT>
ModelInputHandler::ModelInputHandler(
    afw::image::MaskedImage<PixelT> const & image, afw::geom::Point2D const & center, 
    afw::detection::Footprint const & region, int growFootprint,
    afw::image::MaskPixel badPixelMask, bool usePixelWeights,
    double maxBadPixelFraction
) {
    if (growFootprint) {
        _footprint = afw::detection::growFootprint(region, growFootprint);
    } else {
        _footprint = boost::make_shared<afw::detection::Footprint>(region);
    }
    double originalArea = _footprint->getArea();
    _footprint->intersectMask(*image.getMask(), badPixelMask);
    if ((1.0 - _footprint->getArea() / originalArea) > maxBadPixelFraction) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "Too many masked pixels."
        );
    } 
    if (_footprint->getArea() <= 0) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "Model fit contains no usable pixels."
        );
    }
    _data = ndarray::allocate(_footprint->getArea());
    _weights = ndarray::allocate(_footprint->getArea());
    afw::detection::flattenArray(*_footprint, image.getImage()->getArray(), _data, image.getXY0());
    afw::detection::flattenArray(*_footprint, image.getVariance()->getArray(), _weights, image.getXY0());
    if (!usePixelWeights) {
        _weights.asEigen().setConstant(_weights.asEigen().mean());
    }
    _weights.asEigen<Eigen::ArrayXpr>().operator=(_weights.asEigen<Eigen::ArrayXpr>().sqrt().inverse());
    _data.asEigen<Eigen::ArrayXpr>() *= _weights.asEigen<Eigen::ArrayXpr>();
    initCoords(_x, _y, *_footprint, center);
}

template <typename PixelT>
ModelInputHandler::ModelInputHandler(
    afw::image::MaskedImage<PixelT> const & image, afw::geom::Point2D const & center,
    std::vector<afw::geom::ellipses::Ellipse> const & ellipses,
    afw::detection::Footprint const & region, int growFootprint,
    afw::image::MaskPixel badPixelMask, bool usePixelWeights,
    double maxBadPixelFraction
) {
    if (growFootprint) {
        _footprint = afw::detection::growFootprint(region, growFootprint);
    } else {
        _footprint = boost::make_shared<afw::detection::Footprint>(region);
    }
    _footprint = mergeFootprintWithEllipses(*_footprint, ellipses, image.getBBox(afw::image::PARENT));
    double originalArea = _footprint->getArea();
    _footprint->intersectMask(*image.getMask(), badPixelMask);
    if ((1.0 - _footprint->getArea() / originalArea) > maxBadPixelFraction) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "Too many masked pixels."
        );
    } 
    if (_footprint->getArea() <= 0) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "Model fit contains no usable pixels."
        );
    }
    _data = ndarray::allocate(_footprint->getArea());
    _weights = ndarray::allocate(_footprint->getArea());
    afw::detection::flattenArray(*_footprint, image.getImage()->getArray(), _data, image.getXY0());
    afw::detection::flattenArray(*_footprint, image.getVariance()->getArray(), _weights, image.getXY0());
    if (!usePixelWeights) {
        _weights.asEigen().setConstant(_weights.asEigen().mean());
    }
    _weights.asEigen<Eigen::ArrayXpr>().operator=(_weights.asEigen<Eigen::ArrayXpr>().sqrt().inverse());
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
        afw::image::Image<T> const & image, afw::geom::Point2D const & center, \
        std::vector<afw::geom::ellipses::Ellipse> const & ellipses,     \
        afw::detection::Footprint const & region, int growFootprint);   \
    template ModelInputHandler::ModelInputHandler(                      \
        afw::image::MaskedImage<T> const & image, afw::geom::Point2D const & center, \
        afw::geom::Box2I const & region, afw::image::MaskPixel badPixelMask, \
        bool usePixelWeights, double maxBadPixelFraction);              \
    template ModelInputHandler::ModelInputHandler(                      \
        afw::image::MaskedImage<T> const & image, afw::geom::Point2D const & center, \
        afw::detection::Footprint const & region, int growFootprint,    \
        afw::image::MaskPixel badPixelMask, bool usePixelWeights, double maxBadPixelFraction); \
    template ModelInputHandler::ModelInputHandler(                      \
        afw::image::MaskedImage<T> const & image, afw::geom::Point2D const & center, \
        std::vector<afw::geom::ellipses::Ellipse> const & ellipses,     \
        afw::detection::Footprint const & region, int growFootprint,    \
        afw::image::MaskPixel badPixelMask, bool usePixelWeights, double maxBadPixelFraction)

INSTANTIATE(float);
INSTANTIATE(double);

}}}} // namespace lsst::meas::extensions::multiShapelet
