import RiverObs.RiverReachWriter

def write(reach_collection, fout_node, fout_reach, driver='ESRI Shapefile'):
    """
    Writes shapefiles
    """
    reach_variables = list(reach_collection[0].metadata.keys())

    # get node output variables from populated attributes of river_reaches
    node_variables = list(reach_collection[0].__dict__.keys())
    node_variables.remove('ds')
    node_variables.remove('metadata')

    # Write shapefiles
    writer = RiverObs.RiverReachWriter(
        reach_collection, node_variables, reach_variables)

    writer.write_nodes_ogr(fout_node, driver)
    writer.write_reaches_ogr(fout_reach, driver)

