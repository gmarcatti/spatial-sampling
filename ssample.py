# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 13:26:04 2021

@author: Adm
"""

##---- Importar bibliotecas necessárias
import os
from osgeo import ogr
from shapely.wkt import loads
from shapely.geometry import MultiPolygon, MultiPoint, Polygon, Point, asPoint
import numpy as np
from sklearn.cluster import KMeans
from scipy.spatial import Voronoi

##---- Caminho para pasta de trabalho
path = r'C:\Dados\spatial_sampling'

##---- Arquivos shapefile de entrada
talhao = r'\talhao_pol.shp'
amostra = r'\amostra_pts.shp'

##############################################################################
#####--------------------------------------------------------------------#####
#####------ Funções input - output --------------------------------------#####
#####--------------------------------------------------------------------#####
##############################################################################

##--------------------------------------------------------------------------##
##---- Importar feição espacial
##--------------------------------------------------------------------------##

def read_shp(in_shape, m_id):
    in_driver = ogr.GetDriverByName("ESRI Shapefile")
    in_data = in_driver.Open(in_shape, 0)
    in_layer = in_data.GetLayer()
    spatialRef = in_layer.GetSpatialRef()
    mID_shy = []
    in_list = []
    for n in range(in_layer.GetFeatureCount()):
        feature = in_layer.GetFeature(n)
        mID_shy.append(feature.GetField(m_id))
        wkt_feature = loads(feature.geometry().ExportToWkt())
        in_list.append(wkt_feature)
    return((in_list, mID_shy, spatialRef))


##--------------------------------------------------------------------------##
##---- Função exportar feição espacial
##--------------------------------------------------------------------------##

def write_shp(shy, out_shape, id_shy, geom_type, SR = None):
    # geom_type=ogr.wkbPoint ou geom_type=ogr.wkbPolygon
    out_driver = ogr.GetDriverByName("ESRI Shapefile")
    out_data = out_driver.CreateDataSource(out_shape)
    out_layer = out_data.CreateLayer("shy", geom_type=geom_type)
    out_layer.CreateField(ogr.FieldDefn("m_id", ogr.OFTInteger))
    if id_shy is None: id_shy = range(1, len(shy)+1)
    for i, p in enumerate(shy):
        feature = ogr.Feature(out_layer.GetLayerDefn())
        feature.SetField("m_id", id_shy[i])
        if geom_type == 1:
            x, y = p.coords.xy
            wkt = "POINT(%f %f)" %  (float(x[0]) , float(y[0]))
            shape = ogr.CreateGeometryFromWkt(wkt)
        elif geom_type == 3:
            x, y = p.exterior.coords.xy
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for pts in zip(x,y):
                x_i, y_i = pts
                ring.AddPoint(x_i, y_i)
            shape = ogr.Geometry(ogr.wkbPolygon)
            shape.AddGeometry(ring)
        else:
            raise ValueError("Geometria deve ser ogr.wkbPoint ou ogr.wkbPolygon")
        feature.SetGeometry(shape)
        out_layer.CreateFeature(feature)
        feature = None
    out_data = None
    if not SR is None:
        base = os.path.basename(out_shape)
        nome = os.path.splitext(base)[0] + '.prj'
        path = os.path.dirname(out_shape)
        file = open(path + '\\' + nome, 'w')
        file.write(SR.ExportToWkt())
        file.close()


# in_pol = path + talhao
# pol_list, id_pol, SR = read_shp(in_pol, "cod_talhao")
# write_shp(pol_list,  path + r'\talhao_pol_out.shp', range(len(id_pol)), ogr.wkbPolygon, SR)

# in_pts = path + amostra
# pts_list, id_pts, SR  = read_shp(in_pts, "id_amostra")
# write_shp(pts_list,  path + r'\id_amostra_out.shp', range(len(id_pts)), ogr.wkbPoint, SR)

##############################################################################
#####--------------------------------------------------------------------#####
#####------ Funções esquemas amostrais ----------------------------------#####
#####--------------------------------------------------------------------#####
##############################################################################


##--------------------------------------------------------------------------##
##---- Gerar grid amostral - amostras regulares
##--------------------------------------------------------------------------##

def amostra_regula(pol_shy, num_points):
    def grid_regular(pol_shy, num_points):
        num_step = np.int(np.sqrt(pol_shy.area) / np.sqrt(num_points))
        min_x, min_y, max_x, max_y = pol_shy.bounds
        min_x, min_y = (min_x-num_step, min_y-num_step)
        max_x, max_y = (max_x+num_step, max_y+num_step)
        offset = np.random.uniform(0, num_step)
        x = (np.arange(min_x, max_x, step = num_step)) + offset
        y = (np.arange(min_y, max_y, step = num_step)) + offset
        xx,yy = np.meshgrid(x,y)
        xx = xx.reshape((np.prod(xx.shape),))
        yy = yy.reshape((np.prod(yy.shape),))
        #- Converter coordenadas array numpy para pontos shapely
        #points_list = [Point(i) for i in zip(xx,yy)]
        points_list = list(map(asPoint, zip(xx, yy)))
        points_inter = [pts for pts in points_list if pts.within(pol_shy)] 
        #points_shy = MultiPoint(points_inter)
        return(points_inter)
        
    points_shy = grid_regular(pol_shy, num_points)
    while len(points_shy) != num_points:
        points_shy = grid_regular(pol_shy, num_points)
        #print(len(points_shy))
    
    return points_shy

in_pol = path + talhao
pol_list, id_pol, SR = read_shp(in_pol, "cod_talhao")
pol_shy = MultiPolygon(pol_list) # pol_list[0] #

inten = 5 # uma amostra para cada 5 (cinco) hectares
num_points = np.int(pol_shy.area / 10000 / inten)
pontos_regular = amostra_regula(pol_shy, num_points)

MultiPoint(pontos_regular)


##--------------------------------------------------------------------------##
##---- Gerar grid amostral - amostras hexagonais com k-médias
##--------------------------------------------------------------------------##

def amostra_hex(pol, num_points, space_point = 10):
    ##---- Função para converter poligono para pontos
    def pol2pts(pol, space_point):
        min_x, min_y, max_x, max_y = pol.bounds
        x_pol = np.arange(min_x, max_x, step = space_point)
        y_pol = np.arange(min_y, max_y, step = space_point)
        xx_pol, yy_pol = np.meshgrid(x_pol, y_pol)
        xx_pol = xx_pol.reshape((np.prod(xx_pol.shape),))
        yy_pol = yy_pol.reshape((np.prod(yy_pol.shape),))
        pts_list = []
        # pixelWidth = num_step
        for pts in zip(xx_pol, yy_pol):    
            if Point(pts).within(pol):
                pts_list.append(Point(pts))
        return(pts_list)
    
    pts_shy = pol2pts(pol, space_point = space_point)
    X = [[p.x, p.y] for p in pts_shy] 
    kmeans = KMeans(n_clusters = num_points, random_state=0).fit(X)
    pts_coords = kmeans.cluster_centers_
    points_list = list(map(asPoint, pts_coords))
    return points_list


# in_pol = path + talhao
# pol_list, id_pol, SR = read_shp(in_pol, "cod_talhao")
# pol_shy = MultiPolygon(pol_list) # pol_list[0] #

# inten = 5 # uma amostra para cada 5 (cinco) hectares
# num_points = np.int(pol_shy.area / 10000 / inten)
# pontos_hexagono = amostra_hex(pol_shy, num_points)

# MultiPoint(pontos_hexagono)


##############################################################################
#####--------------------------------------------------------------------#####
#####------ Funções auxiliares ------------------------------------------#####
#####--------------------------------------------------------------------#####
##############################################################################

##--------------------------------------------------------------------------##
##---- Gerar os polígonos de voronoi a partir de coordenadas x e y
##--------------------------------------------------------------------------##

def coords2vor(x, y, pol):
    pts_shy = [Point(p) for p in zip(x, y)]
    coords_arr = np.array(list(zip(x, y)))
    min_x, min_y, max_x, max_y = pol.bounds
    max_dist = Point(min_x, min_y).distance(Point(max_x, max_y)) * 3
    min_x, min_y = (min_x-max_dist, min_y-max_dist)
    max_x, max_y = (max_x+max_dist, max_y+max_dist)
    bounds = [[max_x,max_y], [min_x,max_y], [max_x,min_y], [min_x,min_y]]
    coords_arr = np.append(coords_arr, bounds, axis = 0)
    vor = Voronoi(coords_arr)
    pol = pol.buffer(0)
    vor_list = []
    for region in vor.regions:
        if not -1 in region:
            polygon_vor = Polygon([vor.vertices[i] for i in region])
            if any([True for p in pts_shy if p.within(polygon_vor)]):
                vor_list.append(polygon_vor)
    vor_shy = [vor_i.intersection(pol) for vor_i in vor_list]
    return(vor_shy)


# XY = [[p.x, p.y] for p in pontos_regular]
# xx, yy = zip(*XY)
# vor_shy = coords2vor(xx, yy, pol_shy)
# print(np.array([pol.area for pol in vor_shy]) / 10000)
# MultiPolygon(vor_shy)

# XY = [[p.x, p.y] for p in pontos_hexagono]
# xx, yy = zip(*XY)
# vor_shy = coords2vor(xx, yy, pol_shy)
# print(np.array([pol.area for pol in vor_shy]) / 10000)
# MultiPolygon(vor_shy)


##############################################################################
#####--------------------------------------------------------------------#####
#####------ Executar procedimentos --------------------------------------#####
#####--------------------------------------------------------------------#####
##############################################################################

in_pol = path + talhao
pol_list, id_pol, SR = read_shp(in_pol, "cod_talhao")
#pol_shy = MultiPolygon(pol_list) # pol_list[0] #
projeto_vor = []
projeto_pts = []
inten = 1 # uma amostra para cada 5 (cinco) hectares
for i, pol_i in enumerate(pol_list):
    #pol_i = pol_list[0]
    num_points = np.int(pol_i.area / 10000 / inten)
    pontos_hexagono = amostra_hex(pol_i, num_points)
    XY = [[p.x, p.y] for p in pontos_hexagono]
    xx, yy = zip(*XY)
    vor_shy = coords2vor(xx, yy, pol_i)
    projeto_vor.extend(vor_shy)
    projeto_pts.extend(pontos_hexagono)
    vor_area = np.array([pol.area for pol in vor_shy]) / 10000
    print(id_pol[i], "- N:", len(vor_area), "; Média:", round(vor_area.mean(), 2), "; Desvio padrão:", round(vor_area.std(), 4))


# Salver resultados em shapefiles
pts_out = r'\pts_hex1.shp'
pts_out_shape = path + pts_out
write_shp(projeto_pts, pts_out_shape, None, ogr.wkbPoint, SR)

pol_out = r'\voronoi_hex1.shp'
pol_out_shape = path + pol_out
write_shp(projeto_vor, pol_out_shape, None, ogr.wkbPolygon, SR)
