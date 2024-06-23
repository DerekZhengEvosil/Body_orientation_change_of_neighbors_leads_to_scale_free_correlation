from scipy.spatial.distance import pdist, squareform
from get_ellipse_visual_area import get_ellipse_visual_area
from pylab import *


class swarmAlgorithm():
    def __init__(self, r_agent, align_type, select_type, tau, wall_agent, wall_agent_vel, cycTime,
                 fish_BL, w, dist_thresh, config):
        self.align_type = align_type
        self.dist_thresh = dist_thresh
        self.tau = tau
        self.select_type = select_type
        self.wall_agent = wall_agent
        self.wall_agent_vel = wall_agent_vel
        self.fish_BL = fish_BL
        self.w = w
        self.wall_diameter = 0.05
        self.r_agent = r_agent
        self.cycTime = cycTime
        self.v0 = 15
        self.v0_min = 2
        self.rigid_rep_fov = np.pi
        self.maxRotRate = 25 * 1.91
        self.maxRotRate_low = 25 * 1.91

        self.k_align = config.getfloat('Simulation', 'k_align')
        self.k_rep = config.getfloat('Simulation', 'k_rep')
        self.Drep = config.getfloat('Simulation', 'Drep')
        self.dis_agent_agent_to_rep = config.getfloat('Simulation', 'dis_agent_agent_to_rep')
        self.k_rigid_agent_rep = config.getfloat('Simulation', 'k_rigid_agent_rep')
        self.k_r_sense = np.inf
        self.r_sense = self.k_r_sense * self.r_agent

    def calculateVision(self, robotsList, pos_BL, v):
        heading = np.arctan2(v[1, :], v[0, :])
        angular_area, adjacency_matrix, visual_field, polygon_param_cell, polygon_param_belonging_cell, tangent_pt_obj_cart, metric_distance_center = \
            get_ellipse_visual_area(self.w, len(robotsList), self.fish_BL / self.fish_BL, heading, pos_BL, l=1.,
                                    threshold=0., vis_thresh=0.0, dist_thresh=self.dist_thresh,
                                    non_overlap_angles=None)
        adjacency_dist_thresh = {}
        for robot_id in range(1, len(robotsList) + 1):
            belong_tmp = polygon_param_belonging_cell[robot_id]
            neighbor_id_list = [id_tmp for id_tmp in belong_tmp.values()]
            if any(neighbor_id_list):
                adjacency_dist_thresh[robot_id] = np.unique(neighbor_id_list)
            else:
                adjacency_dist_thresh[robot_id] = []

        return angular_area, adjacency_matrix, visual_field, polygon_param_cell, \
               polygon_param_belonging_cell, tangent_pt_obj_cart, metric_distance_center, adjacency_dist_thresh

    def calculateProjection(self, robotsList, angular_area, polygon_param_cell, polygon_param_belonging_cell):
        ellipse_vis_param = {}
        vis_angle = {}
        if np.any(angular_area):
            for r in robotsList:
                polygon_param = polygon_param_cell[r.id + 1]
                polygon_param_belonging = polygon_param_belonging_cell[r.id + 1]
                tangent_point_dis = np.zeros(len(robotsList))
                vis_angle_tmp = zeros(len(robotsList))
                belonging_id_tmp = {}
                polygon_info_tmp = {}
                for pp in range(len(polygon_param)):
                    pp = pp + 1  # index starts from 1
                    belonging_id = polygon_param_belonging[pp]
                    belonging_id_tmp[pp] = belonging_id
                    polygon_param_tmp = polygon_param[pp]
                    polygon_param_p1 = polygon_param_tmp[0]
                    polygon_param_p2 = polygon_param_tmp[1]
                    polygon_param_focal = polygon_param_tmp[2]
                    polygon_info_tmp[pp] = [polygon_param_p1, polygon_param_p2, polygon_param_focal]
                    focal_2_p1 = polygon_param_p1 - polygon_param_focal
                    focal_2_p1 = focal_2_p1 / np.linalg.norm(focal_2_p1)
                    focal_2_p2 = polygon_param_p2 - polygon_param_focal
                    focal_2_p2 = focal_2_p2 / np.linalg.norm(focal_2_p2)
                    tangent_point_dis[belonging_id] = tangent_point_dis[belonging_id] + \
                                                      np.linalg.norm(polygon_param_p1 - polygon_param_p2)
                    dot_temp = np.dot(focal_2_p1, focal_2_p2)
                    if dot_temp > 1.0:
                        dot_temp = 1.0
                    elif dot_temp < -1.0:
                        dot_temp = -1.0
                    vis_angle_tmp[belonging_id] = vis_angle_tmp[belonging_id] + math.acos(dot_temp)
                ellipse_vis_param[r.id] = tangent_point_dis
                vis_angle[r.id] = vis_angle_tmp
        else:
            ellipse_vis_param = {}
            vis_angle = {}
        return ellipse_vis_param, vis_angle

    def calculateAction(self, robotsList, pos, v, all_past_posvel, angular_area, past_vision, adjacency_dist_thresh, informed_id):
        heading = np.arctan2(v[1, :], v[0, :])
        dist_xy = squareform(pdist(pos.T))
        align_mate = []
        old_heading = heading
        if self.align_type == "selective":
            if self.select_type == "BOC":
                for r in robotsList:
                    neighbor = adjacency_dist_thresh[r.id + 1]
                    if any(neighbor):
                        Nneig = neighbor.shape[0]
                        projection_rate = np.zeros(len(robotsList))
                        # accumulate the projection rate within the period
                        for qq in range(1, len(past_vision)):
                            projection_rate = projection_rate + (
                                        abs(past_vision[qq][r.id] - past_vision[qq - 1][r.id]) / self.cycTime)
                        cij = np.ones(Nneig)
                        iDir = v[:, r.id]
                        rij = pos[:, neighbor] - pos[:, r.id].reshape(2, 1)
                        dist_ij = np.linalg.norm(rij, axis=0)
                        nij = rij / dist_ij
                        front_preference = np.power(((1 + np.dot(iDir, nij)) / 2), 0)
                        cij = projection_rate[neighbor] * front_preference
                        salience = np.vstack((neighbor[np.argsort(cij)[::-1]], cij[np.argsort(cij)[::-1]])).T
                        salience[:, 1] = salience[:, 1] / sum(salience[:, 1])
                        align_mate.append(salience)
                    else:
                        align_mate.append(np.array([r.id]))
            elif self.select_type == "random":
                for r in robotsList:
                    neighbor = adjacency_dist_thresh[r.id + 1]
                    if len(neighbor) == 0:
                        align_mate.append(np.array([[r.id, 1]]))
                    else:
                        random_select = np.random.randint(0, len(neighbor))
                        align_mate.append(np.array([[neighbor[random_select], 1]]))

        desSpeed = self.v0 * np.ones([1, len(robotsList)])
        for r in robotsList:
            dist_agent = dist_xy[r.id, :]
            dist_agent[r.id] = np.nan
            dist_agent_wall = np.sqrt(np.sum(np.square(self.wall_agent - pos[:, r.id].reshape((2, 1))), axis=0))
            iDir = v[:, r.id]
            rij = pos[:, :] - pos[:, r.id].reshape(2, 1)
            dist_ij = np.linalg.norm(rij, axis=0)
            dist_ij[r.id] = nan
            nij = rij / dist_ij
            relative_heading = np.arccos(np.dot(iDir, nij))
            see_neighbor = \
                np.where((relative_heading > -self.rigid_rep_fov / 2) & (relative_heading < self.rigid_rep_fov / 2))[0]

            if self.align_type == "selective":
                if len(align_mate[r.id]) > 0:
                    # adapt selection
                    # index2 = align_mate[r.id][:, 1]
                    # index3 = np.where(np.round(cumsum(index2), 6) >= self.adap_percent)[0][0]
                    # index1 = align_mate[r.id][0:index3 + 1, 0]
                    # index2 = index2[0:index3 + 1]
                    # max selection
                    if len(align_mate[r.id].shape) == 2:
                        index2 = np.max(align_mate[r.id][:, 1])
                        index1 = align_mate[r.id][0, 0]
                    else:
                        index2 = 1
                        index1 = np.array([r.id])
                else:
                    index1 = []
                    index2 = []
                neighbor_and_self = v[:, np.append(r.id, index1.astype(np.int))]
                aij = np.append(1, 1)
                # aij[1:] = aij[1:] / sum(aij[1:])
                part_vel = np.nansum(aij * neighbor_and_self, axis=1)
                part_vel = part_vel / np.linalg.norm(part_vel)

            rep_neighbor_agent = np.where(dist_agent < self.Drep)[0]
            if np.size(rep_neighbor_agent):
                nij = (pos[:, r.id].reshape((2, 1)) - pos[:, rep_neighbor_agent]) / dist_agent[rep_neighbor_agent]
                part_pos = (self.Drep - dist_agent[rep_neighbor_agent]) * nij
                part_pos = np.sum(part_pos, axis=1)
                part_pos = part_pos / np.linalg.norm(part_pos)
            else:
                part_pos = np.array([0, 0])

            rigidrep_neighbor_agent = np.intersect1d(np.where(dist_agent < self.dis_agent_agent_to_rep)[0],
                                                     see_neighbor)
            if np.size(rigidrep_neighbor_agent):
                nij = (pos[:, r.id].reshape((2, 1)) - pos[:, rigidrep_neighbor_agent]) / dist_agent[
                    rigidrep_neighbor_agent]
                part_rigid_rep = (self.Drep - dist_agent[rigidrep_neighbor_agent]) * nij
                part_rigid_rep = np.sum(part_rigid_rep, axis=1)
                part_rigid_rep = part_rigid_rep / np.linalg.norm(part_rigid_rep)
            else:
                part_rigid_rep = np.array([0, 0])

            part = np.array([0, 0])
            if np.size(rigidrep_neighbor_agent):
                part = self.k_rigid_agent_rep * part_rigid_rep
                desSpeed[0, r.id] = self.v0_min
            else:
                part = self.k_align * part_vel + self.k_rep * part_pos

            part = part / np.linalg.norm(part)
            desDir = part
            curDir = v[:, r.id]
            rotateDir = np.arcsin(curDir[0] * desDir[1] - curDir[1] * desDir[0])  # x1y2 – x2y1
            dot_temp = desDir[0] * curDir[0] + desDir[1] * curDir[1]
            if dot_temp > 1.0:
                dot_temp = 1.0
            elif dot_temp < -1.0:
                dot_temp = -1.0
            rotateAnlge = np.arccos(dot_temp) * 180 / np.pi
            rotRate = rotateAnlge / self.cycTime
            multiple_TR = 0
            if desSpeed[0, r.id] == self.v0:
                if rotRate > self.maxRotRate:
                    rotRate = self.maxRotRate
                    multiple_TR = math.floor(rotRate / 1.91)
                else:
                    rotRate = rotRate
                    multiple_TR = math.floor(rotRate / 1.91)
            elif desSpeed[0, r.id] == self.v0_min:
                if rotRate > self.maxRotRate:
                    rotRate = self.maxRotRate
                    multiple_TR = math.floor(rotRate / 1.91)
                else:
                    rotRate = rotRate
                    multiple_TR = math.floor(rotRate / 1.91)

            if rotateDir > 0:
                multiple_TR = multiple_TR
            else:
                multiple_TR = -1 * multiple_TR

            robotLeftspeed = desSpeed[0, r.id] + multiple_TR
            robotRightspeed = desSpeed[0, r.id] - multiple_TR
            r.desSpeed = desSpeed[0, r.id]
            r.rotateRate = rotRate
            r.multiple_TR = multiple_TR
            r.leftSpeed = robotLeftspeed
            r.rightSpeed = robotRightspeed

        return robotsList
