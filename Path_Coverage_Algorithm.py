import math
import os
import sys
from enum import IntEnum
import matplotlib.pyplot as plt
import numpy as np

class Grid:

    def __init__(self, width, height, resolution,center_x, center_y, init_val=0.0):
        self.width = width
        self.height = height
        self.resolution = resolution
        self.center_x = center_x
        self.center_y = center_y

        self.left_lower_x = self.center_x - self.width / 2.0 * self.resolution
        self.left_lower_y = self.center_y - self.height / 2.0 * self.resolution

        self.ndata = self.width * self.height
        self.data = [init_val] * self.ndata


    def  ConvertIndex_xy_position_xy(self, x_pos, y_pos):
        x_ind = self. CalculateIndexvalue_position(x_pos, self.left_lower_x, self.width)
        y_ind = self. CalculateIndexvalue_position(y_pos, self.left_lower_y, self.height)

        return x_ind, y_ind

    def ValueFromIndex(self, x_ind, y_ind):
        grid_ind = self.GridIndex_from_CordIndex(x_ind, y_ind)

        if 0 <= grid_ind < self.ndata:
            return self.data[grid_ind]
        else:
            return None



    def set_value_from_xy_pos(self, x_pos, y_pos, val):

        x_ind, y_ind = self. ConvertIndex_xy_position_xy(x_pos, y_pos)

        if (not x_ind) or (not y_ind):
            return False  # NG

        flag = self.SettingValues_Index(x_ind, y_ind, val)

        return flag

    def SettingValues_Index(self, x_ind, y_ind, val):
        if (x_ind is None) or (y_ind is None):
            return False, False

        grid_ind = int(y_ind * self.width + x_ind)

        if 0 <= grid_ind < self.ndata:
            self.data[grid_ind] = val
            return True  # OK
        else:
            return False  # NG

    def PolygonValues(self, pol_x, pol_y, val, inside=True):
        # making ring polygon
        if (pol_x[0] != pol_x[-1]) or (pol_y[0] != pol_y[-1]):
            pol_x.append(pol_x[0])
            pol_y.append(pol_y[0])

        # setting value for all grid
        for x_ind in range(self.width):
            for y_ind in range(self.height):
                x_pos, y_pos = self.CentralGridIndexValues(x_ind, y_ind)
                flag = self.check_inside_polygon(x_pos, y_pos, pol_x, pol_y)

                if flag is inside:
                    self.SettingValues_Index(x_ind, y_ind, val)

    def GridIndex_from_CordIndex(self, x_ind, y_ind):
        grid_ind = int(y_ind * self.width + x_ind)
        return grid_ind

    def CentralGridIndexValues(self, x_ind, y_ind):
        x_pos = self.calc_grid_central_xy_position_from_index(x_ind, self.left_lower_x)
        y_pos = self.calc_grid_central_xy_position_from_index(y_ind, self.left_lower_y)

        return x_pos, y_pos

    def calc_grid_central_xy_position_from_index(self, index, lower_pos):
        return lower_pos + index * self.resolution + self.resolution / 2.0

    def  CalculateIndexvalue_position(self, pos, lower_pos, max_index):
        ind = int(np.floor((pos - lower_pos) / self.resolution))
        if 0 <= ind <= max_index:
            return ind
        else:
            return None

    def CheckingOccupation(self, xind, yind, occupied_val=1.0):

        val = self.ValueFromIndex(xind, yind)

        if val is None or val >= occupied_val:
            return True
        else:
            return False

    def GridExpand(self):
        xinds, yinds = [], []

        for ix in range(self.width):
            for iy in range(self.height):
                if self.CheckingOccupation(ix, iy):
                    xinds.append(ix)
                    yinds.append(iy)

        for (ix, iy) in zip(xinds, yinds):
            self.SettingValues_Index(ix + 1, iy, val=1.0)
            self.SettingValues_Index(ix, iy + 1, val=1.0)
            self.SettingValues_Index(ix + 1, iy + 1, val=1.0)
            self.SettingValues_Index(ix - 1, iy, val=1.0)
            self.SettingValues_Index(ix, iy - 1, val=1.0)
            self.SettingValues_Index(ix - 1, iy - 1, val=1.0)

    @staticmethod
    def check_inside_polygon(iox, ioy, x, y):

        npoint = len(x) - 1
        inside = False
        for i1 in range(npoint):
            i2 = (i1 + 1) % (npoint + 1)

            if x[i1] >= x[i2]:
                min_x, max_x = x[i2], x[i1]
            else:
                min_x, max_x = x[i1], x[i2]
            if not min_x < iox < max_x:
                continue

            tmp1 = (y[i2] - y[i1]) / (x[i2] - x[i1])
            if (y[i1] + tmp1 * (iox - x[i1]) - ioy) > 0.0:
                inside = not inside

        return inside

    def PlotGrid(self, ax=None):

        grid_data = np.reshape(np.array(self.data), (self.height, self.width))
        if not ax:
            fig, ax = plt.subplots()
        heat_map = ax.pcolor(grid_data, cmap="BuPu", vmin=0.0, vmax=1.0)
        plt.axis("equal")
        plt.show()

        return heat_map

do_animation = True


class SweepSearcher:
    class SweepDirection(IntEnum):
        UP = 1
        DOWN = -1

    class MovingDirection(IntEnum):
        RIGHT = 1
        LEFT = -1

    def __init__(self,
                 moving_direction, sweep_direction, x_inds_goal_y, goal_y):
        self.moving_direction = moving_direction
        self.sweep_direction = sweep_direction
        self.turing_window = []
        self.update_turning_window()
        self.xinds_goaly = x_inds_goal_y
        self.goaly = goal_y

    def movingGrid(self, cxind, cyind, gmap):
        nxind = self.moving_direction + cxind
        nyind = cyind

        # found safe grid
        if not gmap.CheckingOccupation(nxind, nyind,
                                                 occupied_val=0.5):
            return nxind, nyind
        else:  # occupied
            ncxind, ncyind = self.SaferGrid(cxind, cyind, gmap)
            if (ncxind is None) and (ncyind is None):
                # moving backward
                ncxind = -self.moving_direction + cxind
                ncyind = cyind
                if gmap.CheckingOccupation(ncxind, ncyind):
                    # moved backward, but the grid is occupied by obstacle
                    return None, None
            else:
                # keep moving until end
                while not gmap.CheckingOccupation(ncxind + self.moving_direction, ncyind,occupied_val=0.5):
                    ncxind += self.moving_direction
                self.SwapMove()
            return ncxind, ncyind

    def SaferGrid(self, cxind, cyind, gmap):

        for (d_x_ind, d_y_ind) in self.turing_window:

            next_x_ind = d_x_ind + cxind
            next_y_ind = d_y_ind + cyind

            # found safe grid
            if not gmap.CheckingOccupation(next_x_ind, next_y_ind,
                                                     occupied_val=0.5):
                return next_x_ind, next_y_ind

        return None, None

    def SearchDone(self, grid):
        for ix in self.xinds_goaly:
            if not grid.CheckingOccupation(ix, self.goaly,
                                                         occupied_val=0.5):
                return False

        # all lower grid is occupied
        return True

    def update_turning_window(self):
        # turning window definition
        # robot can move grid based on it.
        self.turing_window = [
            (self.moving_direction, 0.0),
            (self.moving_direction, self.sweep_direction),
            (0, self.sweep_direction),
            (-self.moving_direction, self.sweep_direction),
        ]

    def SwapMove(self):
        self.moving_direction *= -1
        self.update_turning_window()

    def StartGrid(self, grid):
        x_inds = []
        y_ind = 0
        if self.sweep_direction == self.SweepDirection.DOWN:
            x_inds, y_ind =SearchEdges(
                grid, from_upper=True)
        elif self.sweep_direction == self.SweepDirection.UP:
            x_inds, y_ind =SearchEdges(
                grid, from_upper=False)

        if self.moving_direction == self.MovingDirection.RIGHT:
            return min(x_inds), y_ind
        elif self.moving_direction == self.MovingDirection.LEFT:
            return max(x_inds), y_ind

        raise ValueError("self.moving direction is invalid ")


def StartPosition_Direction(ox, oy):
    # find sweep_direction
    max_dist = 0.0
    vec = [0.0, 0.0]
    sweep_start_pos = [0.0, 0.0]
    for i in range(len(ox) - 1):
        dx = ox[i + 1] - ox[i]
        dy = oy[i + 1] - oy[i]
        d = np.hypot(dx, dy)

        if d > max_dist:
            max_dist = d
            vec = [dx, dy]
            sweep_start_pos = [ox[i], oy[i]]

    return vec, sweep_start_pos


def convert_grid_coordinate(ox, oy, sweep_vec, sweep_start_posi):
    tx = [ix - sweep_start_posi[0] for ix in ox]
    ty = [iy - sweep_start_posi[1] for iy in oy]

    th = math.atan2(sweep_vec[1], sweep_vec[0])

    c = np.cos(-th)
    s = np.sin(-th)

    rx = [ix * c - iy * s for (ix, iy) in zip(tx, ty)]
    ry = [ix * s + iy * c for (ix, iy) in zip(tx, ty)]

    return rx, ry


def convert_global_coordinate(x, y, sweep_vec, sweep_start_posi):
    th = math.atan2(sweep_vec[1], sweep_vec[0])
    c = np.cos(th)
    s = np.sin(th)

    tx = [ix * c - iy * s for (ix, iy) in zip(x, y)]
    ty = [ix * s + iy * c for (ix, iy) in zip(x, y)]

    rx = [ix + sweep_start_posi[0] for ix in tx]
    ry = [iy + sweep_start_posi[1] for iy in ty]

    return rx, ry


def SearchEdges(grid, from_upper=False):
    yind = None
    xinds = []

    if from_upper:
        xrange = range(grid.height)[::-1]
        yrange = range(grid.width)[::-1]
    else:
        xrange = range(grid.height)
        yrange = range(grid.width)

    for iy in xrange:
        for ix in yrange:
            if not grid.CheckingOccupation(ix, iy):
                yind = iy
                xinds.append(ix)
        if yind:
            break
    # print(xinds, yind)
    return xinds, yind


def SetupGrid(ox, oy, resolution, sweep_direction, offset_grid=10):
    width = math.ceil((max(ox) - min(ox)) / resolution) + offset_grid
    height = math.ceil((max(oy) - min(oy)) / resolution) + offset_grid
    center_x = (np.max(ox)+np.min(ox))/2.0
    center_y = (np.max(oy)+np.min(oy))/2.0

    grid = Grid(width, height, resolution, center_x, center_y)
    grid.PolygonValues(ox, oy, 1.0, inside=False)
    grid.GridExpand()

    x_inds_goal_y = []
    goal_y = 0
    if sweep_direction == SweepSearcher.SweepDirection.UP:
        x_inds_goal_y, goal_y =SearchEdges(grid,from_upper=True)
    elif sweep_direction == SweepSearcher.SweepDirection.DOWN:
        x_inds_goal_y, goal_y =SearchEdges(grid,from_upper=False)

    return grid, x_inds_goal_y, goal_y


def SweepPath(sweep_searcher, grid, grid_search_animation=False):
    # search start grid
    cxind, cyind = sweep_searcher.StartGrid(grid)
    if not grid.SettingValues_Index(cxind, cyind, 0.5):
        print("Cannot find start grid")
        return [], []

    x, y = grid.CentralGridIndexValues(cxind, cyind)
    px, py = [x], [y]

    fig, ax = None, None
    if grid_search_animation:
        fig, ax = plt.subplots()
        # for stopping simulation with the esc key.
        fig.canvas.mpl_connect('key_release_event',
                               lambda event: [
                                   exit(0) if event.key == 'escape' else None])

    while True:
        cxind, cyind = sweep_searcher.movingGrid(cxind, cyind, grid)

        if sweep_searcher.SearchDone(grid) or (cxind is None or cyind is None):
            print("Done")
            break
        x, y = grid.CentralGridIndexValues(cxind, cyind)
        px.append(x)
        py.append(y)

        grid.SettingValues_Index(cxind, cyind, 0.5)
        if grid_search_animation:
            grid.PlotGrid(ax=ax)
            plt.pause(1.0)
    grid.PlotGrid()

    return px, py


def PathPlanning(ox, oy, resolution,moving_direction=SweepSearcher.MovingDirection.RIGHT,sweeping_direction=SweepSearcher.SweepDirection.UP,):
    sweep_vec, sweep_start_posi = StartPosition_Direction(ox, oy)
    rox, roy = convert_grid_coordinate(ox, oy, sweep_vec, sweep_start_posi)
    gmap, xinds_goaly, goaly = SetupGrid(rox, roy, resolution,sweeping_direction)
    sweep_searcher = SweepSearcher(moving_direction, sweeping_direction,xinds_goaly, goaly)
    px, py = SweepPath(sweep_searcher, gmap)
    rx, ry = convert_global_coordinate(px, py, sweep_vec, sweep_start_posi)
    print("Path length:", len(rx))
    return rx, ry


def PathAnimation(ox, oy, resolution):
    px, py = PathPlanning(ox, oy, resolution)
    if do_animation:
        for ipx, ipy in zip(px, py):
            plt.cla()
            plt.gcf().canvas.mpl_connect('key_release_event',lambda event: [exit(0) if event.key == 'escape' else None])
            plt.plot(ox, oy, "b")
            plt.plot(px, py, "black")
            plt.plot(ipx, ipy, "or")
            plt.axis("equal")
            plt.grid(True)
            plt.pause(0.01)
    plt.cla()
    plt.plot(ox, oy, "b")
    plt.plot(px, py, "black")
    plt.axis("equal")
    plt.pause(10)
    plt.close()


def main():
    print("start!!")
    
    ox = [20.0, 20.0,35.0,55.0,68.0,90.0,104.0,108.0,100.0,78.0,55.0,20.0]
    oy = [20.0, 50.0,80.0,105.0,109.0,109.0,98.0,75.0,40.0,33.0,20.0,20.0]
    resolution = 2.5
    PathAnimation(ox, oy, resolution)

    # ox = [25.0,30.0,38.0,28.0,28.0,45.0,52.0,65.0,70.0,79.0,85.0,96.0,86.0,90.0,80.0,96.0,100.0,105.0,108.0,109.0,106.0,102.0,95.0,78.0,66.0,60.0,55.0,46.0,40.0,38.0,25.0]
    # oy = [60.0,65.0,78.0,85.0,100,104.0,96.0,90.0,92.0,100.0,112.0,105.0,92.0,80.0,76.0,66.0,80.0,83.0,71.0,60,56.0,40.0,38.0,49.0,51.0,40.0,38.0,56.0,58.0,56.0,60.0]
    # print(len(ox))
    # print(len(oy))
    # resolution = 1
    # PathAnimation(ox, oy, resolution)

#     ox = [0.0, 20.0, 50.0, 200.0, 130.0, 40.0, 0.0]
#     oy = [0.0, -80.0, 0.0, 30.0, 60.0, 80.0, 0.0]
#     resolution = 3.4
#     PathAnimation(ox, oy, resolution)

    plt.show()
    print("Solution Found!!!")


if __name__ == '__main__':
    main()
