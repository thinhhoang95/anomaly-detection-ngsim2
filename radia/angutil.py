def roundang(ang): # -180 to 180
    if ang <= -180:
            ang = ang + 360
    if ang >= 180:
        ang = ang - 360
    return ang

def roundang2(angs): # same as roundang function but takes 2 parameters
    new_angs = ()
    for ang in angs:
        if ang <= -180:
                ang = ang + 360
        if ang >= 180:
            ang = ang - 360
        new_angs = new_angs + (ang,)
    return new_angs


def conv_std_to_360(ang): # from 0 to 360
    if ang>=0 and ang<=180:
        ang = ang
    if ang<0 and ang>-180:
        ang = ang + 360
    return ang

def degdiff(ang1, ang2): # take difference and then perform roundang function -> always < 180deg
    diff = ang1 - ang2
    if diff <= -180:
            diff = diff + 360
    if diff >= 180:
        diff = diff - 360
    return diff

def sort_increase(angs): # first param will be reached first while traversing counter clockwise on the complex unit circle
    ang1=angs[0]
    ang2=angs[1]
    candidate_value = -1000
    if degdiff(ang1, ang2) > candidate_value:
        max_ang = ang1
        min_ang = ang2
        candidate_value = degdiff(ang1, ang2)
    if degdiff(ang2, ang1) > candidate_value:
        max_ang = ang2
        min_ang = ang1
        candidate_value = degdiff(ang2, ang1)
    return (min_ang, max_ang)