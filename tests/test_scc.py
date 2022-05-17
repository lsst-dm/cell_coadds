# Test runner
# NOT a unit test

# import skyproj
import healsparse as hsp
import healpy
import sys
sys.path.append("/software/lsstsw/stack_20220125/conda/miniconda3-py38_4.9.2/envs/lsst-scipipe/lib/python3.8/site-packages/")
import numpy as np

import lsst.daf.butler as dafButler
import lsst.geom as geom
import lsst.sphgeom as sphgeom
from lsst.pipe.tasks.coaddBase import makeSkyInfo
from lsst.pipe.tasks.makeSkyMap import MakeSkyMapTask, MakeSkyMapConfig

import lsst.skymap
import lsst.cell_coadds

# butler = dafButler.Butler("/repo/main_20210215/", collections=["HSC/runs/RC2/w_2021_50/DM-32987"])
butler = dafButler.Butler("/repo/dc2/", collections=["2.2i/runs/test-med-1/w_2021_48/DM-32707"])

where_string = "skymap='DC2_cells_v1' AND tract=3828 AND patch=45 AND band='r'"
where_dataId = butler.registry.expandDataId({"skymap": "DC2_cells_v1", "tract":3828, "patch":45, "band": 'r'})
calexps = [butler.getDirectDeferred(ref) for ref in butler.registry.queryDatasets("calexp", dataId=where_dataId).expanded()]
calexps.sort(key=lambda handle: handle.dataId)

skymap0 = butler.get("skyMap", dataId={"skymap": "DC2_cells_v1"})
skyInfo0 = makeSkyInfo(skymap0, tractId=3828, patchId=45) # (1,4)

config = lsst.cell_coadds.CoaddInCellsConfig()
task = lsst.cell_coadds.CoaddInCellsTask(config=config)

multipleCellCoadd = task.run({"calexps":calexps}, skyInfo0, where_dataId)
for idx, scc in enumerate(multipleCellCoadd.cells):
    mi = scc.inner.asMaskedImage()
    mi.writeFits(f"/scratch/kannawad/cellCoadd_{idx}.fits")