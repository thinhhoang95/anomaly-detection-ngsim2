import numpy as np
from radia import angutil as ang
from radia import radar as rd
import math
import sys, io
import logging
from radia import v2x as v2x

def isNeighborRayOff(pos1, ang1, pos2):
    # The given vehicle will be neighbor to this vehicle if it is located
    # in side the radar cone ahead (1) and surrounding radar sphere (2)
    # (2) surrounding radar sphere
    pos1x = np.array(list(pos1))
    pos2x = np.array(list(pos2))
    posx = pos2x - pos1x

    dist = np.linalg.norm(posx)
    
    if dist<rd.SRADAR_RADIUS:
        sphere_radar_visible = True
    else:
        sphere_radar_visible = False
    
    if dist<rd.FRADAR_RADIUS:
        front_radar_visible = True
        # bearing = np.arctan2(posx[1],posx[0])/np.pi*180 # global bearing 
        # logging.log(25, 'Bearing: ' + str(bearing))
        # delta_alpha = ang.degdiff(bearing, ang1)
        # logging.log(25, 'Delta Alpha: ' + str(delta_alpha))
        # if np.abs(delta_alpha) < rd.FRADAR_HALF_ANGLE:
        #     front_radar_visible = True
        # else:
        #     front_radar_visible = False
    else:
        front_radar_visible = False

    return (front_radar_visible, sphere_radar_visible)

def isAcquaintance(pos1, pos2):
    # The given vehicle will be acquaintance to this vehicle if it is located
    # in side the communicable circle

    pos1x = np.array(list(pos1))
    pos2x = np.array(list(pos2))
    posx = pos2x - pos1x

    dist = np.linalg.norm(posx)
    
    if dist<v2x.COMM_RANGE:
        acquaintance = True
    else:
        acquaintance = False
    return acquaintance

def establishComm(host, globs):
    pos1 = host.posx
    for obj in globs: 
        if obj.id == host.id:
            continue
        pos2 = obj.posx
        if isAcquaintance(pos1, pos2):
            host.addAcquaintance(obj.id)
        

def getCorners(center, bearing):
    xC = center[0]
    yC = center[1]
    bearing = bearing * np.pi / 180 # convert degrees to radians
    d = math.sqrt(rd.VEH_LENGTH**2 + rd.VEH_WIDTH**2)
    alpha = np.arctan2(rd.VEH_WIDTH, rd.VEH_LENGTH)
    xA = xC + 0.5 * d * np.cos(bearing + alpha)
    yA = yC + 0.5 * d * np.sin(bearing + alpha)
    xB = xC + 0.5 * d * np.cos(bearing - alpha)
    yB = yC + 0.5 * d * np.sin(bearing - alpha)
    xD = xC - 0.5 * d * np.cos(bearing + alpha)
    yD = yC - 0.5 * d * np.sin(bearing + alpha)
    xE = xC - 0.5 * d * np.cos(bearing - alpha)
    yE = yC - 0.5 * d * np.sin(bearing - alpha)
    return ((xA,yA), (xB, yB), (xD, yD), (xE, yE))

def getObjectRays(corners):
    ang1 = np.arctan2(corners[0][1],corners[0][0]) / np.pi * 180
    ang2 = np.arctan2(corners[1][1],corners[1][0]) / np.pi * 180
    ang3 = np.arctan2(corners[2][1],corners[2][0]) / np.pi * 180
    ang4 = np.arctan2(corners[3][1],corners[3][0]) / np.pi * 180
    candidate_value = -1000
    if ang.degdiff(ang1, ang2) > candidate_value:
        max_ang = ang1
        min_ang = ang2
        candidate_value = ang.degdiff(ang1, ang2)
    if ang.degdiff(ang2, ang1) > candidate_value:
        max_ang = ang2
        min_ang = ang1
        candidate_value = ang.degdiff(ang2, ang1)
    if ang.degdiff(ang1, ang3) > candidate_value:
        max_ang = ang1
        min_ang = ang3
        candidate_value = ang.degdiff(ang1, ang3)
    if ang.degdiff(ang3, ang1) > candidate_value:
        max_ang = ang3
        min_ang = ang1
        candidate_value = ang.degdiff(ang3, ang1)
    if ang.degdiff(ang1, ang4) > candidate_value:
        max_ang = ang1
        min_ang = ang4
        candidate_value = ang.degdiff(ang1, ang4)
    if ang.degdiff(ang4, ang1) > candidate_value:
        max_ang = ang4
        min_ang = ang1
        candidate_value = ang.degdiff(ang4, ang1)
    if ang.degdiff(ang2, ang3) > candidate_value:
        max_ang = ang2
        min_ang = ang3
        candidate_value = ang.degdiff(ang2, ang3)
    if ang.degdiff(ang3, ang2) > candidate_value:
        max_ang = ang3
        min_ang = ang2
        candidate_value = ang.degdiff(ang3, ang2)
    if ang.degdiff(ang2, ang4) > candidate_value:
        max_ang = ang2
        min_ang = ang4
        candidate_value = ang.degdiff(ang2, ang4)
    if ang.degdiff(ang4, ang2) > candidate_value:
        max_ang = ang4
        min_ang = ang2
        candidate_value = ang.degdiff(ang4, ang2)
    if ang.degdiff(ang3, ang4) > candidate_value:
        max_ang = ang3
        min_ang = ang4
        candidate_value = ang.degdiff(ang3, ang4)
    if ang.degdiff(ang4, ang3) > candidate_value:
        max_ang = ang4
        min_ang = ang3
        candidate_value = ang.degdiff(ang4, ang3)
        
    return (min_ang, max_ang)


# Experimental: ray tracing must imply that the object must significantly smaller than the radar angle

def frontRadarRayTrace(host, globs): # globs: global objects
    host.neighbors.clear()
    pos1 = host.posx
    ang1 = host.orientation
    # Sort the global objects with respect to the Euclidean distance
    sorted_globs = sorted(globs, key=lambda x: np.linalg.norm(x.posx - pos1))
    # Perform ray tracing for each individual objects in the global array of vehicles
    # All will be in host's coordinate system
    ray_min = ang.roundang(- rd.FRADAR_FULL_ANGLE/2) + 90
    ray_max = ang.roundang(rd.FRADAR_FULL_ANGLE/2) + 90
    ray_min, ray_max = ang.sort_increase((ray_min, ray_max))
    free_angles = [(ray_min, ray_max)]
    for obj in sorted_globs: 
        if (obj.id == host.id):
            continue # skip the current iteration
        center = obj.posx 
        # Eliminate too far away targets for speeded up calculations
        front, _ = isNeighborRayOff(pos1, ang1, center) 
        if host.id == 679 and obj.id == 690:
            print('Host: ', host.posx, ', Object: ', obj.posx)
        if not front:
            continue

        # bearing = obj.orientation
        relative_center = obj.posx - host.posx
        relative_bearing = ang.roundang(obj.orientation - host.orientation)
        # Find 4 corners of the objects
        corners = getCorners(relative_center, relative_bearing) # in the host's coordinate system
        objectRays = ang.roundang2(getObjectRays(corners))
        objectRays = ang.sort_increase(objectRays)
        new_free_angles = []
        for free_angle in free_angles: 
            free_angle = ang.sort_increase(free_angle)
            if ang.degdiff(objectRays[0],free_angle[0])>=0 and ang.degdiff(objectRays[1],free_angle[1])<=0:
            # if objectRays[0]>=free_angle[0] and objectRays[1]<=free_angle[1]:
                # The object is observable, repartition the free_angles array
                new_free_angles.append((free_angle[0], objectRays[0]))
                new_free_angles.append((objectRays[1], free_angle[1]))
                host.addNeighbor(obj.id)
            elif ang.degdiff(objectRays[1],free_angle[1])>=0 and ang.degdiff(objectRays[0],free_angle[0])>=0 and ang.degdiff(objectRays[0],free_angle[1])<=0:
            # elif objectRays[1]>=free_angle[0] and objectRays[0]<=free_angle[0]:
                # The object is partially present on the left
                #print('Object presents on the left')
                host.addNeighbor(obj.id) # maybe we should add a threshold on 
                new_free_angles.append((free_angle[0], objectRays[0]))
            elif ang.degdiff(objectRays[0],free_angle[0])<=0 and ang.degdiff(objectRays[1],free_angle[1]) and ang.degdiff(objectRays[1],free_angle[0])>=0:
            # elif objectRays[1]>=free_angle[1] and objectRays[0]<=free_angle[1]:
                # The object is partially present on the right
                #print('Object presents on the right')
                host.addNeighbor(obj.id)
                new_free_angles.append((objectRays[1], free_angle[1]))
            elif ang.degdiff(objectRays[0],free_angle[1])>=0 or ang.degdiff(objectRays[1],free_angle[0])<=0:
            # elif objectRays[0]>=free_angle[1] or objectRays[1]<=free_angle[0]:
                # The object is outside of the region of interest
                new_free_angles.append((free_angle[0], free_angle[1]))
            else:
                sys.exit('Error: Unknown state while performing ray tracing') 
        free_angles = new_free_angles
    return (host, free_angles)

def surroundRadarRayTrace(host, globs): # similar to frontRadarRayTrace
    pos1 = host.posx
     # Sort the global objects with respect to the Euclidean distance
    sorted_globs = sorted(globs, key=lambda x: np.linalg.norm(x.posx - pos1))
    free_angles = [(0, 179.999), (-179.999, 0)]
    
    for obj in sorted_globs:
        if (obj.id == host.id):
            continue
        center = obj.posx
        _, sphere = isNeighborRayOff(pos1, 0, center)
        if not sphere:
            continue
        bearing = obj.orientation
        corners = getCorners(center, bearing)
        objectRays = ang.roundang2(getObjectRays(corners))
        objectRays = ang.sort_increase(objectRays)
        new_free_angles = []
        last_free_angle = ()
        adjacent_index = 0
        
        for free_angle in free_angles:
            if adjacent_index > 0:
                adjacent_index = adjacent_index - 1
            free_angle = ang.sort_increase(free_angle)
            if ang.degdiff(objectRays[0],free_angle[0])>=0 and ang.degdiff(objectRays[1],free_angle[1])<=0:
            # if objectRays[0]>=free_angle[0] and objectRays[1]<=free_angle[1]:
                # The object is observable, repartition the free_angles array
                new_free_angles.append((free_angle[0], objectRays[0]))
                new_free_angles.append((objectRays[1], free_angle[1]))
                host.addNeighbor(obj.id)
            elif ang.degdiff(objectRays[1],free_angle[0])>=0 and ang.degdiff(objectRays[0],free_angle[0])<=0:
            # elif objectRays[1]>=free_angle[0] and objectRays[0]<=free_angle[0]:
                # The object is partially present on the right
                new_free_angles.append((objectRays[1], free_angle[1]))
                # Special case: if object simultaneously presents on partial left and right, it is due to its location at 0 and 180 deg
                if len(last_free_angle)>0 and adjacent_index == 1 and ang.degdiff(free_angle[0], last_free_angle[1])>-1:
                    host.addNeighbor(obj.id)
                else:
                    adjacent_index = adjacent_index + 2 
            elif ang.degdiff(objectRays[1],free_angle[1])>=0 and ang.degdiff(objectRays[0],free_angle[1])<=0:
            # elif objectRays[1]>=free_angle[1] and objectRays[0]<=free_angle[1]:
                # The object is partially present on the left
                new_free_angles.append((free_angle[0], objectRays[0]))
                # Special case: if object simultaneously presents on partial left and right, it is due to its location at 0 and 180 deg
                if len(last_free_angle)>0 and adjacent_index == 1 and ang.degdiff(free_angle[0], last_free_angle[1])>-1:
                    host.addNeighbor(obj.id)
                else:
                    adjacent_index = adjacent_index + 2 
            elif ang.degdiff(objectRays[0],free_angle[1])>=0 or ang.degdiff(objectRays[1],free_angle[0])<=0:
            # elif objectRays[0]>=free_angle[1] or objectRays[1]<=free_angle[0]:
                # The object is outside of the region of interest
                new_free_angles.append((free_angle[0], free_angle[1]))
            else:
                sys.exit('Error: Unknown state while performing ray tracing') 
            last_free_angle = free_angle
        free_angles = new_free_angles
    return (host, free_angles)
