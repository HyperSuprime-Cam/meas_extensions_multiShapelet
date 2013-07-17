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

#include <list>
#include <algorithm>

#include "boost/format.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/meas/extensions/multiShapelet/MultiGaussianRegistry.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

namespace {

typedef std::pair<std::string,MultiGaussian> RegistryItem;
typedef std::list<RegistryItem> RegistryList;

RegistryList & getRegistryList() {
    static RegistryList it;
    return it;
}

struct CompareRegistryItem {

    bool operator()(RegistryItem const & item) const { return item.first == name; }

    explicit CompareRegistryItem(std::string const & name_) : name(name_) {}

    std::string const & name;
};

} // anonymous

MultiGaussian const & MultiGaussianRegistry::lookup(std::string const & name) {
    RegistryList & l = getRegistryList();
    RegistryList::iterator i = std::find_if(l.begin(), l.end(), CompareRegistryItem(name));
    if (i == l.end()) {
        throw LSST_EXCEPT(
            pex::exceptions::NotFoundException,
            (boost::format("MultiGaussian with name '%s' not found in registry.") % name).str()
        );
    }
    MultiGaussian const & result = i->second;
    if (i != l.begin()) {
        l.splice(l.begin(), l, i);
    }
    return result;
}

void MultiGaussianRegistry::insert(std::string const & name, MultiGaussian const & multiGaussian) {
    RegistryList & l = getRegistryList();
    RegistryList::iterator i = std::find_if(l.begin(), l.end(), CompareRegistryItem(name));
    if (i != l.end()) {
        i->second = multiGaussian;
    } else {
        l.push_back(std::make_pair(name, multiGaussian));
    }
}

void MultiGaussianRegistry::insert(
    std::string const & name,
    ndarray::Array<double const,1> const & fluxes,
    ndarray::Array<double const,1> const & radii,
    bool normalize
) {
    if (fluxes.getSize<0>() != radii.getSize<0>()) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("flux array size (%d) does not match radius array size (%d)") %
             fluxes.getSize<0>() % radii.getSize<0>()).str()
        );
    }
    MultiGaussian multiGaussian;
    double totalFlux = 0.0;
    for (int n = 0; n < fluxes.getSize<0>(); ++n) {
        multiGaussian.add(GaussianComponent(fluxes[n], radii[n]));
        totalFlux += fluxes[n];
    }
    if (normalize) {
        for (int n = 0; n < fluxes.getSize<0>(); ++n) {
            multiGaussian[n].flux /= totalFlux;
        }
    }
    insert(name, multiGaussian);
}

}}}} // namespace lsst::meas::extensions::multiShapelet
