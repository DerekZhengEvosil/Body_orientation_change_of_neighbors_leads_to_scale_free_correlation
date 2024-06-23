import numpy as np
import itertools as it
import sys
# import matplotlib.pyplot as plt
# from matplotlib.patches import Ellipse, Polygon, FancyArrowPatch

def _get_tangent_point_parameter(w, r, theta, phi, fish_BL):
    '''calculates where the tangent points lie on the ellipse, return the corresponding angles,
    these can be translated in to coordinates via using the function
    ellipse_point_from_parameter()
    '''
    w = w / 2.0
    main_axis = fish_BL/2.
    aa = np.sqrt(-2.0 * main_axis * main_axis * w * w + (main_axis * main_axis + w * w) * r * r +
                 (w * w - main_axis * main_axis) * r * r * np.cos(2.0 * (theta - phi))) / np.sqrt(2.0)  # 公式中的gamma

    bb = w * r * np.cos(theta - phi) - main_axis * w  # 公式中的beta

    psi1 = 2.0 * np.arctan2(aa - main_axis * r * np.sin(theta - phi), bb)
    psi2 = -2.0 * np.arctan2(aa + main_axis * r * np.sin(theta - phi), bb)
    return [psi1, psi2]


def _ellipse_point_from_parameter(r, theta, phi, psi, w, fish_BL):
    # calculates cartesian coordinates for a point on an ellipse
    # with long axis 1, short axis w, ellipse center at r,theta
    # that is given by the ellipse parameter psi
    l = fish_BL/2
    x = r * np.cos(theta) + l * np.cos(phi) * np.cos(psi) + w * l * np.sin(phi) * np.sin(psi)
    y = r * np.sin(theta) + l * np.sin(phi) * np.cos(psi) - w * l * np.cos(phi) * np.sin(psi)
    return [x, y]



def _line_intersect(x1a, y1a, x1b, y1b, x2a, y2a, x2b, y2b):
    # finds the point where two lines intersect.
    # lines are given by two points each:
    # line 1 goes through (x1a,y1a) and (x1b, y1b)
    # the same goes for line2

    point = None
    if (x1b - x1a) != 0:
        if (x2b - x2a) != 0:
            m1 = (y1b - y1a) / (x1b - x1a)
            m2 = (y2b - y2a) / (x2b - x2a)
            if m1 != m2:
                b1 = y1a - m1 * x1a
                b2 = y2a - m2 * x2a
                x = (b2 - b1) / (m1 - m2)
                y = m1 * x + b1
                point = [x, y]
            else:
                print('lines are parallel')
        else:
            x = x2b
            m1 = (y1b - y1a) / (x1b - x1a)
            b1 = y1a - m1 * x1a
            y = m1 * x + b1
            point = [x, y]
    else:
        if (x2b - x2a) != 0:
            x = x1b
            m2 = (y2b - y2a) / (x2b - x2a)
            b2 = y2a - m2 * x2a
            y = m2 * x + b2
            point = [x, y]
        else:
            print('lines are parallel')
    if point != None:
        point = np.array(point)
    return point


def _cast_to_pm_pi(a):
    '''Casts any (radian) angle to the
            equivalent in the interval (-pi, pi)'''
    b = (a + np.pi) % (2. * np.pi)
    b -= np.pi
    return b


def _get_closest_id(r, out, n):
    """ used to find the closest intersection point on a ray emitted from and ellipses eye,
        r is numpy array with indices jklm as follows:
        j: which intersection [2],
        k: on which ellipse [n],
        l: for which viewer [n],
        m: for which ray [2(n-1)]"""
    out_r = np.zeros(shape=out.shape, dtype=float)
    for j, k in it.product(range(n), range((n - 1) * 2)):
        if np.isnan(r[:, :, j, k]).all():
            out[j, k] = np.nan
            out_r[j, k] = np.nan
        else:
            out[j, k] = np.nanargmin(r[:, :, j, k], axis=1)[1]
            out_r[j, k] = np.nanmin(r[:, :, j, k])
    return out, out_r


def _get_ellipse_line_intersection_points(eyes, tps, w, fish_BL):
    ''' given two points of the line (eyes and tp) calculates
        the points at which this line intersects with an ellipse
        of length 1 and width w with center at the origin and
        orientation along the positive x-axis,
        returns points as 2x2 array,
        index1: x/y,
        index2: which intersection point,
        if only 1 intersections found both entries are equal,
        if no intersections are found, entries are np.nan'''
    x1 = eyes[0]
    y1 = eyes[1]
    x2 = tps[0]
    y2 = tps[1]
    a = fish_BL/2.
    b = w / 2.
    dd = ((x2 - x1) ** 2 / (a ** 2) + (y2 - y1) ** 2 / (b ** 2))
    ee = (2. * x1 * (x2 - x1) / (a ** 2) + 2. * y1 * (y2 - y1) / (b ** 2))
    ff = (x1 ** 2 / (a ** 2) + y1 ** 2 / (b ** 2) - 1.)
    determinant = ee ** 2 - 4. * dd * ff
    float_epsilon = 0.00001
    zeromask = abs(determinant) >= 1000. * float_epsilon
    determinant *= zeromask
    determinant[np.where(determinant < 0)] = np.nan
    t = (np.array([(-ee - np.sqrt(determinant)) / (2. * dd),
                   (-ee + np.sqrt(determinant)) / (2. * dd)]))
    mask = np.array(t > 0., dtype=float)
    mask[mask == 0.] = np.nan
    x = mask * (x1 + (x2 - x1) * t)
    y = mask * (y1 + (y2 - y1) * t)
    return np.array([x, y])


def _remove_self_intersections(inters, n):
    ''' used to remove intersections of ray emitted from ellipse i's eye and intersecting with
        ellipse i's boundary when detecting all intersections of those rays with all other ellipses,
        inters is array of interception points with indices ijklm
        i: x/y [2],
        j: which intersection [2],
        k: on which ellipse [n],
        l: for which viewer [n],
        m: for which ray [2(n-1)]'''

    for i in range(n):
        inters[:, :, i, i, :] = np.nan
    return inters



def _subjpol_to_objcart(r, theta, pos, phi):
    # takes in a point, r, theta from the polar coordinates
    # with center at pos and orientation phi
    # returns a point in cartesian coordinates (same
    # coordinate system that pos is given in)
    # phi = _cast_to_pm_pi(phi)
    rot_mat_back = np.array([[np.cos(-phi), np.sin(-phi)], [-np.sin(-phi), np.cos(-phi)]])
    pt_subj = [r * np.cos(theta), r * np.sin(theta)]
    pt_obj = np.dot(rot_mat_back, pt_subj) + pos
    return pt_obj


    
def get_ellipse_visual_area(w, n, fish_BL, phi, ori_pos, l=1., threshold=0., vis_thresh=0.0, dist_thresh=np.inf, non_overlap_angles=None):
    n = int(n)
    tp_subj_pol = []
    tp_obj_cart = []
    pos_center = ori_pos
    # metric_distance_center = np.zeros([n, n])
    z_center = np.array([[complex(p[0], p[1]) for p in pos_center.T]])
    metric_distance_center = abs(z_center.T - z_center)
    pos = ori_pos - np.array([-fish_BL * l / 2.0 * np.cos(phi) , -fish_BL * l / 2.0 * np.sin(phi)])
    # rename some variables for convenience
    phi_m = np.array([phi, ] * n).transpose()
    x = pos[0]
    y = pos[1]
    x_center = pos_center[0]
    y_center = pos_center[1]

    # calculate the relative positions of i to j in coordinate
    # system with origin in the eye of j
    rel_x = x_center.reshape(len(x_center), 1) - x  # entry(ij)=pos(i)-pos(j)
    rel_y = y_center.reshape(len(y_center), 1) - y
    theta = np.arctan2(rel_y, rel_x)
    z = np.array([[complex(p[0], p[1]) for p in pos.T]])
    z_center = np.array([[complex(p[0], p[1]) for p in pos_center.T]])
    r = abs(z_center.T - z)
    # indices ij: abs(z_center(i)-z(j)), j is observer, i target

    # to avoid errors in further calc. result for these will be set manually
    np.fill_diagonal(metric_distance_center, float('NaN'))
    np.fill_diagonal(r, float("NaN"))

    psi = _get_tangent_point_parameter(w, r, theta, phi_m, fish_BL)

    for p in psi:
        # calculate tangent point from psi in local polar coordinates
        pt_subj_pol = _ellipse_point_from_parameter(r, theta, phi_m, p, w, fish_BL)
        # print("p_shape{}".format(p[~np.isnan(p)].shape))
        # print("r_shape{}".format(r[~np.isnan(r)].shape))
        z_pt_subj_pol = pt_subj_pol[0] + 1j * pt_subj_pol[1]
        # print("pt_subj_pol[0]_shape{}".format(pt_subj_pol[0][~np.isnan(pt_subj_pol[0])].shape))
        # print("pt_subj_pol[1]_shape{}".format(pt_subj_pol[1][~np.isnan(pt_subj_pol[1])].shape))
        nan_check = np.arctan2(pt_subj_pol[1], pt_subj_pol[0])
        # print("nan_check_shape{}".format(nan_check[~np.isnan(nan_check)].shape))
        theta_tp = _cast_to_pm_pi(np.arctan2(pt_subj_pol[1], pt_subj_pol[0]) - phi)
        r_tp = abs(z_pt_subj_pol)
        np.fill_diagonal(r_tp, 0.0)
        # print("theta_tp_shape{}".format(theta_tp[~np.isnan(theta_tp)].shape))
        tp_subj_pol.append(np.array([r_tp, theta_tp]))  # focal的眼睛相对于这个切点的距离和角度
        # transform tangent points to cartesian global coordinates
        pt_obj_cart = pt_subj_pol + np.array([np.array([pos[0], ] * n), np.array( \
            [pos[1], ] * n)])
        np.fill_diagonal(pt_obj_cart[0], 0.0)
        np.fill_diagonal(pt_obj_cart[1], 0.0)
        tp_obj_cart.append(pt_obj_cart)
    tangent_pt_subj_pol = np.array(tp_subj_pol)
    tangent_pt_obj_cart = np.array(tp_obj_cart)

    # get ray angles for each ellipse
    angles = tangent_pt_subj_pol[:, 1].flatten(order='f')
    nan_angle = angles[~np.isnan(angles)]
    if nan_angle.size != 2 * (n - 1) * n:
        return np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
    # print("angles_shape{}".format(angles[~np.isnan(angles)].shape))
    angles = np.sort(nan_angle.reshape(2 * (n - 1), n, order='f').T)
    assert np.logical_and(angles.all() <= np.pi, angles.all() >= -np.pi), 'angles are not in pm pi interval'
    between_angles = _cast_to_pm_pi(
        np.diff(angles, append=(2. * np.pi + angles[:, 0]).reshape(n, 1), axis=1) / 2. + angles)

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #  transformation of angles for the calculation of intersection points
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # transform the local angles into points in global cartesian coordinates
    phi_hlp = np.repeat(phi.reshape(n, 1), 2 * (n - 1), axis=1)
    transf_betw_ang = between_angles + phi_hlp
    raypoints = np.array([np.cos(transf_betw_ang), np.sin(transf_betw_ang)]) + np.tile(pos, (
        (n - 1) * 2, 1, 1)).transpose(1, 2, 0)

    # here we need to transform the raypoints from global coordinates to local
    # ones of the ellipse that we want to check of intersections
    # (in a manner that will set up a nested for loop)
    raypoints = np.tile(raypoints, (n, 1, 1, 1)).transpose(1, 0, 2, 3)
    # indices: x/y ,N repetitions (in which coordinate system), focalid (seen from which eye), raypoints (which tangent point)

    pos_hlp = np.tile(pos_center, (2 * (n - 1), 1, 1)).transpose(1, 2, 0)
    pos_hlp = np.tile(pos_hlp, (n, 1, 1, 1)).transpose(1, 2, 0, 3)
    # indices: ijkl x/y,id (coordinate syst.=the individual that intersections will be found for), repetition (which eye), repetitions (which tangent point)
    # shifting the raypoints to a coordinate system with origin in the center of the ellipse j (the one that intersections will be found for)
    raypoints -= pos_hlp

    # now go to polar coordinates and rotate the points by -phi,
    # to orient the ellipse j along positive x-axis in the respective
    # coordinate system (this is needed because the function calculating
    # intersections assumes an ellipse at the center with this orientation) #
    r = np.sqrt(raypoints[0] ** 2 + raypoints[1] ** 2)
    theta = np.arctan2(raypoints[1], raypoints[0])
    phi_hlp = np.tile(phi, (n, (n - 1) * 2, 1)).transpose(2, 0, 1)
    theta -= phi_hlp
    # now the transofmration is over
    raypoints = np.array([r * np.cos(theta), r * np.sin(theta)])

    # Now we need to similarly transform the eye positions from
    # global to local (in a manner that will set up a nested for loop)
    # (the id of the viewer ellipse is the second last index, thus
    # the array needs to have repetitions for all other axes)
    eyes = np.tile(pos, (2 * (n - 1), 1, 1)).transpose(1, 2, 0)
    eyes = np.tile(eyes, (n, 1, 1, 1)).transpose(1, 0, 2, 3)
    # shift coordinate system origins
    eyes -= pos_hlp
    # rotate coordinate systems
    r = np.sqrt(eyes[0] ** 2 + eyes[1] ** 2)
    theta = np.arctan2(eyes[1], eyes[0])
    theta -= phi_hlp
    eyes = np.array([r * np.cos(theta), r * np.sin(theta)])
    inters = _get_ellipse_line_intersection_points(eyes, raypoints, w, fish_BL)
    inters = _remove_self_intersections(inters, n)  #
    # indices: [x/y, which intersection, on which ellipse,
    # for which viewer, for which ray]
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++=

    # all intersection points are still in coordinates of
    # the 'on which ellipse' ellipse, transform to global coordinates next:
    # 1. rotate by +phi #
    theta = np.arctan2(inters[1], inters[0]) + phi_hlp
    r = np.sqrt(inters[0] ** 2 + inters[1] ** 2)
    inters = np.array([r * np.cos(theta), r * np.sin(theta)])
    # 2. and shift position of origin
    pos_hlp = np.tile(pos_hlp, (2, 1, 1, 1, 1)).transpose(1, 0, 2, 3, 4)
    inters = inters + pos_hlp

    # in order to decide which intersection point is closest to an
    # ellipse we need to move to the coordinate system of the ellipse
    # which is emitting the rays from its eye (second last index)
    # (we skip the rotation because we are only interested in the
    # distances r anyways)
    pos_hlp = np.tile(pos, (2 * (n - 1), 1, 1)).transpose(1, 2, 0)
    pos_hlp = np.tile(pos_hlp, (n, 1, 1, 1)).transpose(1, 2, 0, 3)
    pos_hlp = np.tile(pos_hlp, (2, 1, 1, 1, 1)).transpose(1, 0, 3, 2, 4)
    # shift to the local coordinates
    inters -= pos_hlp
    # calculate the distances:  # intersection点再次转换到局部坐标中，相对于原点的距离
    r = np.sqrt(inters[0] ** 2 + inters[1] ** 2)

    # Here want to find for each ray emitted from the eye of a viewer ellipse,
    # the id of the closest ellipse it intersects with
    out = np.empty([n, (n - 1) * 2], dtype=float)
    closest_r = []
    closest_id, closest_r = _get_closest_id(r, out, n)
    visual_field = np.stack([closest_id, angles, np.roll(angles, -1, axis=-1)])
    area = np.stack([closest_id, (np.diff(visual_field[1::, :, :], axis=0) % np.pi)[0]])
    angular_area = np.zeros([n, n], dtype=float)
    for i in range(n):
        mask = area[0] == i
        angular_area[i, :] = np.sum(mask * area[1], axis=-1)
    adjacency_matrix = np.array(angular_area > threshold, dtype=int)

    # get the data for drawing results
    segments = visual_field
    tps = tangent_pt_obj_cart
    md = metric_distance_center
    polygon_param_cell = {}
    polygon_param_belonging_cell = {}
    for viewer_id in range(n):
        cnt = 0
        polygon_param = {}
        polygon_param_belonging = {}
        for k in range(2 * (n - 1)):
            if not np.isnan(segments[0, viewer_id, k]):
                i = int(segments[0, viewer_id, k])
                if angular_area[i, viewer_id] > vis_thresh and md[i, viewer_id] < dist_thresh:
                    cnt = cnt + 1
                    hlp_low = _subjpol_to_objcart(md[i, viewer_id], segments[1, viewer_id, k], pos[:, viewer_id],
                                                phi[viewer_id])
                    hlp_high = _subjpol_to_objcart(md[i, viewer_id], segments[2, viewer_id, k], pos[:, viewer_id],
                                                phi[viewer_id])
                    p1 = _line_intersect(hlp_low[0], hlp_low[1], pos[0, viewer_id], pos[1, viewer_id],
                                        tps[0, 0, i, viewer_id], tps[0, 1, i, viewer_id], tps[1, 0, i, viewer_id],
                                        tps[1, 1, i, viewer_id])
                    p2 = _line_intersect(hlp_high[0], hlp_high[1], pos[0, viewer_id], pos[1, viewer_id],
                                        tps[0, 0, i, viewer_id], tps[0, 1, i, viewer_id], tps[1, 0, i, viewer_id],
                                        tps[1, 1, i, viewer_id])
                    
                    polygon_param[cnt] = [p1, p2, pos[:, viewer_id]]
                    polygon_param_belonging[cnt] = i

        polygon_param_cell[viewer_id + 1] = polygon_param
        polygon_param_belonging_cell[viewer_id + 1] = polygon_param_belonging

    # # draw the results
    # fig_show = plt.gcf()
    # plt_ax = plt.gca()
    # zorder = 100
    # viewer_id = 0
    # vis_thresh = 0.0
    # dist_thresh = np.inf
    # alpha = 0.7
    # edgecolor = 'none'
    # edgewidth = 1
    # vf_color = 'darkseagreen'
    # show_eyes = 1
    # eyecolor = 'k'
    # eyesize = 5
    # ellipse_color = "r"
    # ellipses_list = []
    # # draw individual
    # print(pos_center)
    # for i in range(n):
    #     ellipses = Ellipse(pos_center[:, i], 1.0, w,
    #                        _cast_to_pm_pi(phi[i]) * 180.0 / np.pi)
    #     plt_ax.add_artist(ellipses)
    #     ellipses.set_clip_box(plt_ax.bbox)
    #     ellipses.set_facecolor(ellipse_color)
    #     ellipses.set_alpha(alpha)
    #     ellipses.set_edgecolor(edgecolor)
    #     ellipses.set_linewidth(edgewidth)
    #     ellipses.set_zorder(zorder)
    #     ellipses_list.append(ellipses)
    #
    # cnt = 0
    # polygon_param = {}
    # for k in range(2 * (n - 1)):
    #     if not np.isnan(segments[0, viewer_id, k]):
    #         cnt = cnt + 1
    #         i = int(segments[0, viewer_id, k])
    #         # print(i)
    #         if angular_area[i, viewer_id] > vis_thresh and md[i, viewer_id] < dist_thresh:
    #             hlp_low = _subjpol_to_objcart(md[i, viewer_id], segments[1, viewer_id, k], pos[:, viewer_id],
    #                                           phi[viewer_id])
    #             hlp_high = _subjpol_to_objcart(md[i, viewer_id], segments[2, viewer_id, k], pos[:, viewer_id],
    #                                            phi[viewer_id])
    #             p1 = _line_intersect(hlp_low[0], hlp_low[1], pos[0, viewer_id], pos[1, viewer_id],
    #                                  tps[0, 0, i, viewer_id], tps[0, 1, i, viewer_id], tps[1, 0, i, viewer_id],
    #                                  tps[1, 1, i, viewer_id])
    #             p2 = _line_intersect(hlp_high[0], hlp_high[1], pos[0, viewer_id], pos[1, viewer_id],
    #                                  tps[0, 0, i, viewer_id], tps[0, 1, i, viewer_id], tps[1, 0, i, viewer_id],
    #                                  tps[1, 1, i, viewer_id])
    #             visual_area = plt.Polygon([p1, p2, pos[:, viewer_id]])
    #             polygon_param[cnt] = [p1, p2, pos[:, viewer_id]]
    #             plt_ax.add_artist(visual_area)
    #             visual_area.set_facecolor(vf_color)
    #             visual_area.set_alpha(alpha)
    #             ellipses_list.append(visual_area)
    #
    #     plt_ax.set_aspect('equal')
    #     plt_ax.set_xlim(np.amin(pos_center[0]) - 1, np.amax(pos_center[0]) + 1)
    #     plt_ax.set_ylim(np.amin(pos_center[1]) - 1, np.amax(pos_center[1]) + 1)
    #     plt_ax.spines['top'].set_visible(False)
    #     plt_ax.spines['right'].set_visible(False)
    #     plt_ax.tick_params(axis='both', colors='0.5')
    #     plt_ax.spines['bottom'].set_color('0.5')
    #     plt_ax.spines['left'].set_color('0.5')
    #
    # plt.show()

    return angular_area, adjacency_matrix, visual_field, polygon_param_cell, polygon_param_belonging_cell, tangent_pt_obj_cart, metric_distance_center


def _generate_initial_spatial_configuration(state, nn, noise_int, d, kappa, w):
    if state == 'grid':
        n = int(np.floor(np.sqrt(nn)))
        xlen = n
        ylen = n
        number = n * n
        grid_x = np.linspace(d, d * xlen, xlen, endpoint=True)
        grid_y = np.linspace(d, d * ylen, ylen, endpoint=True)
        x, y = np.meshgrid(grid_x, grid_y)
        pos = np.array([x.flatten(), y.flatten()])
        if n < np.sqrt(nn):
            for i in range(nn - number):
                extra = np.array([d * (xlen + 1 + np.floor(i / n)), d * (i % n + 1)]).reshape(2, 1)
                pos = np.hstack([pos, extra])
        orientations = np.random.vonmises(0.0, kappa, nn)
        noise = (np.random.random((2, nn)) - np.ones((2, nn)) * 0.5) * 2.0 * noise_int * d
        pos = pos + noise
        return pos, orientations

    elif state == 'hexagonal':
        d_y = d / np.sqrt(2.)
        n = int(np.floor(np.sqrt(nn)))
        xlen = n
        ylen = n
        number = n * n
        grid_x = np.linspace(d, d * xlen, xlen, endpoint=True)
        grid_y = np.linspace(d_y, d_y * ylen, ylen, endpoint=True)
        x, y = np.meshgrid(grid_x, grid_y)
        x[0:-1:2] += d / 2.
        pos = np.array([x.flatten(), y.flatten()])
        if n < np.sqrt(nn):
            for i in range(nn - number):
                extra = np.array([d * (xlen + 1 + np.floor(i / n)), d_y * (i % n + 1)]).reshape(2, 1)
                pos = np.hstack([pos, extra])
        orientations = np.random.vonmises(0.0, kappa, nn)
        noise_x = (np.random.random((nn)) - np.ones((nn)) * 0.5) * 2.0 * noise_int * d
        noise_y = (np.random.random((nn)) - np.ones((nn)) * 0.5) * 2.0 * noise_int * d_y
        pos[0] += noise_x
        pos[1] += noise_y
        return pos, orientations
    else:
        print("state needs to be either milling or grid or hexagonal")
