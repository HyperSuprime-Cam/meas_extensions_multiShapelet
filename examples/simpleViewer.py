import lsst.meas.extensions.multiShapelet as ms
import lsst.afw.image
import lsst.afw.geom.ellipses
import lsst.shapelet
from  lsst.meas.algorithms.algorithmRegistry import AlgorithmRegistry

CANONICAL_FLUX = "flux.psf"

def makePsfImage(source, bbox):
    """
    Make an image of the multi-shapelet PSF approximation saved in the given source.

    Arguments:
      source ------ An afw.table.SourceRecord with 'multishapelet.psf.*' fields.
      bbox -------- An afw.geom.Box2I giving the location and size of the image to make.
                    Should usually be the same as the bbox of the image returned by calling
                    computeImage on the Psf model at the source's slot-defined centroid.
    
    Returns an afw.image.ImageD.
    """
    ctrl = ms.FitPsfControl()
    model = ms.FitPsfModel(ctrl, source)
    func = model.asMultiShapelet(source.getCentroid())
    image = lsst.afw.image.ImageD(bbox)
    func.evaluate().addToImage(image)
    return image

def removeFluxCorrections(source, algName):
    """
    Given the name of an aperture-corrected and PSF-fit-corrected flux field, return
    uncorrected flux and flux error for comparison to the original pixels.
    """
    flux = source.get(algName + ".flux")
    fluxErr = source.get(algName + ".flux.err")

    # remove aperture correction
    apCorr = source.get("aperturecorrection")
    apCorrErr = source.get("aperturecorrection.err")
    flux1 = flux / apCorr
    fluxErr1 = ((fluxErr**2 - flux**2 * apCorrErr**2) / apCorr**2)**0.5
    
    # remove PSF-fit correction
    canonicalPsfFactor = source.get(CANONICAL_FLUX + ".psffactor") # canonical algorithm on Psf realization
    psfFactor = source.get(algName + ".psffactor") # this algorithm on Psf realization
    scaling = canonicalPsfFactor / psfFactor
    flux2 = flux1 / scaling
    fluxErr2 = fluxErr1 / scaling
    
    return flux2, fluxErr2

def makeGalaxyImages(source, bbox=None):
    """
    Make images of the multi-shapelet models saved in the given source.

    Arguments:
      source ------ An afw.table.SourceRecord with 'multishapelet.psf.*' fields.
      bbox -------- An afw.geom.Box2I giving the location and size of the image to make.
                    If None, a box that includes the 4-r_e ellipse of both profile models
                    will be used.

    Returns three afw.image.ImageDs, ordered (combo, exp, dev), all with the same dimensions.
    """
    haveInputBBox = True
    if bbox is None:
        bbox = lsst.afw.geom.Box2I()
        haveInputBBox = False

    funcs = []

    if source.get("multishapelet.combo.flux.flags"):
        raise ValueError("Source has flags set; bbox may not be sane.")

    psfCtrl = ms.FitPsfControl()
    psfModel = ms.FitPsfModel(psfCtrl, source)
    psfFunc = psfModel.asMultiShapelet()
    psfFunc.normalize()

    for name, config in zip(("exp", "dev"), (ms.FitExponentialConfig, ms.FitDeVaucouleurConfig)):
        ctrl = config().makeControl()
        ctrl.name = "multishapelet." + name
        model = ms.FitProfileModel(ctrl, source)
        model.flux, model.fluxErr = removeFluxCorrections(source, ctrl.name)
        func = model.asMultiShapelet(source.getCentroid())
        func = func.convolve(psfFunc)
        funcs.append(func)

        if not haveInputBBox:
            quadrupole = model.ellipse.convolve(psfModel.ellipse)
            ellipse = lsst.afw.geom.ellipses.Ellipse(quadrupole, source.getCentroid())
            ellipse.getCore().scale(4.0)
            bboxD = ellipse.computeEnvelope()
            bbox.include(lsst.afw.geom.Box2I(bboxD))

    comboFlux, comboFluxErr = removeFluxCorrections(source, "multishapelet.combo")
    devFrac = source.get("multishapelet.combo.devfrac")

    comboImage = lsst.afw.image.ImageD(bbox)
    images = [comboImage]
    for func, weight in zip(funcs, (1.0 - devFrac, devFrac)):
        image = lsst.afw.image.ImageD(bbox)
        func.evaluate().addToImage(image)
        images.append(image)
        func.normalize(weight * comboFlux)
        func.evaluate().addToImage(comboImage)

    return images

