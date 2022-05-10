# THIS SHOULD BE IN DESC REPO USING make_coadd_obs from descwl_coadd

from __future__ import annotations

from typing import Mapping, Tuple

import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
import lsst.sphgeom
import lsst.utils
import numpy as np
from descwl_coadd import make_coadd_obs
#from lsst.cell_coadds import SingleCellCoaddBuilder
from ._cellCoaddBuilder import SingleCellCoaddBuilder, singleCellCoaddBuilderRegistry
from ._image_planes import OwnedImagePlanes
from lsst.pex.config import Field, registerConfigurable
from lsst.pipe.tasks.coaddBase import makeSkyInfo


@registerConfigurable("sccBuilder", singleCellCoaddBuilderRegistry)
class SCCBuilder(SingleCellCoaddBuilder):
    """A concrete class to build single cell coadds"""

    def run(
        self,
        inputs: Mapping[ObservationIdentifiers, Tuple[DeferredDatasetHandle, geom.Box2I]],
        cellInfo: pipeBase.Struct,
    ):
        # import pdb; pdb.set_trace()
        coadd_obs, exp_info = make_coadd_obs(
            exps=[_v[0].get(bbox=_v[1]) for _k, _v in inputs.items()],
            coadd_wcs=cellInfo.wcs,
            coadd_bbox=cellInfo.outer_bbox,
            psf_dims=(41, 41), # self.config.psf_dims,
            rng=np.random.RandomState(12345),
            remove_poisson=True,  # no object poisson noise in sims
        )

        # TODO learn how to save the noise exp as well
        center = coadd_obs.coadd_exp.psf.getAveragePosition()
        image_planes = OwnedImagePlanes(image=coadd_obs.coadd_exp.image,
                                        mask=coadd_obs.coadd_exp.mask,
                                        variance=coadd_obs.coadd_exp.variance,
                                        mask_fractions={})
        return pipeBase.Struct(image_planes=image_planes,
                               psf=coadd_obs.coadd_exp.psf.computeImage(center),
                               inputs=inputs)

    @staticmethod
    def get_noise_exp(exp, rng):
        """Get a noise image based on the input exposure

        TODO gain correct separately in each amplifier, currently
        averaged

        Parameters
        ----------
        exp: `~lsst.afw.image.ExposureF`
            The exposure upon which to base the noise

        Returns
        -------
        noise_exp: `~lsst.afw.image.ExposureF`
            A noise exposure with the same WCS and size as the input exposure.
        """
        signal = exp.image.array
        variance = exp.variance.array.copy()

        use = np.where(np.isfinite(variance) & np.isfinite(signal))

        gains = [amp.getGain() for amp in exp.getDetector().getAmplifiers()]
        mean_gain = np.mean(gains)

        corrected_var = variance[use] - signal[use] / mean_gain

        medvar = np.median(corrected_var)

        noise_image = rng.normal(scale=np.sqrt(medvar), size=signal.shape)

        ny, nx = signal.shape
        nmimage = afwImage.MaskedImageF(width=nx, height=ny)
        assert nmimage.image.array.shape == (ny, nx)

        nmimage.image.array[:, :] = noise_image
        nmimage.variance.array[:, :] = medvar
        nmimage.mask.array[:, :] = exp.mask.array[:, :]

        noise_exp = afwImage.ExposureF(nmimage)
        noise_exp.setPsf(exp.getPsf())
        noise_exp.setWcs(exp.getWcs())
        noise_exp.setFilterLabel(exp.getFilterLabel())
        noise_exp.setDetector(exp.getDetector())

        return noise_exp

    @staticmethod
    def hash_function(seed, tract, patch, band):
        """Generate a hash key given the base seed and metadata"""
        band_map = {"u": 1, "g": 2, "r": 3, "i": 4, "z": 5, "y": 6}
        # Add a linear combination of metadata weighted by prime numbers
        hash_key = seed + 131071 * tract + 524287 * patch + 8388607 * band_map[band]
        return hash_key
