#from shapely.ops import unary_union
import sys
import geopandas as gpd
import logging
import shapely

from ecoshard import taskgraph

logging.basicConfig(
    level=logging.WARNING,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)



AOI_VECTOR_TUPLES = {
    'hybas_sa_lev05_intersect_non-arpa': (
        './data/non-arpa-projected-in-m.gpkg',
        './data/hydrosheds/hybas_sa_lev05_v1c.shp'),
    'hybas_si_lev05_intersect_Arctic_si': (
        './data/Arctic.gpkg',
        './data/hydrosheds/hybas_si_lev05_v1c.shp'),
    'hybas_ar_lev05_intersect_Arctic_ar': (
        './data/Arctic.gpkg',
        './data/hydrosheds/hybas_ar_lev05_v1c.shp'),
    'hybas_sa_lev05_intersect_arpa': (
        './data/arpa-projected-in-m.gpkg',
        './data/hydrosheds/hybas_sa_lev05_v1c.shp'),
    'hybas_sa_lev05_intersect_Colombia': (
        './data/Colombia.gpkg',
        './data/hydrosheds/hybas_sa_lev05_v1c.shp'),
    'hybas_na_lev05_intersect_NGP': (
        './data/NGP.gpkg',
        './data/hydrosheds/hybas_na_lev05_v1c.shp'),
    'hybas_sa_lev05_intersect_Peru': (
        './data/Peru.gpkg',
        './data/hydrosheds/hybas_sa_lev05_v1c.shp'),
    'hybas_na_lev05_intersect_RGRB': (
        './data/RGRB.gpkg',
        './data/hydrosheds/hybas_na_lev05_v1c.shp'),
    'hybas_sa_lev05_intersect_Tapajos': (
        './data/Tapajos.gpkg',
        './data/hydrosheds/hybas_sa_lev05_v1c.shp'),

}


def calculate_intersection(aoi_path, subwatersheds_path, target_path):
    aoi = gpd.read_file(aoi_path)
    subwatersheds = gpd.read_file(subwatersheds_path)

    if aoi.crs != subwatersheds.crs:
        aoi = aoi.to_crs(subwatersheds.crs)

    # Fix possibly invalid geometries by buffering 0
    aoi["geometry"] = aoi["geometry"].buffer(0)
    subwatersheds["geometry"] = subwatersheds["geometry"].buffer(0)

    # Instead of aoi.unary_union (deprecated), use shapely's union_all
    aoi_union = shapely.union_all(aoi["geometry"])

    subwatersheds_intersected = subwatersheds[subwatersheds.intersects(aoi_union)]
    subwatersheds_intersected.to_file(target_path, driver='GPKG')


def main():
    task_graph = taskgraph.TaskGraph('.', len(AOI_VECTOR_TUPLES), 15.0)
    for target_basename, (aoi_path, subwatersheds_path) in AOI_VECTOR_TUPLES.items():
        target_path = f'./data/{target_basename}.gpkg'
        task_graph.add_task(
            func=calculate_intersection,
            args=(aoi_path, subwatersheds_path, target_path),
            target_path_list=[target_path],
            task_name=f'{target_basename}')
    task_graph.join()
    task_graph.close()


if __name__ == '__main__':
    main()
