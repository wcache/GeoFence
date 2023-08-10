from math import radians, cos, sqrt, pi

MAX_LNG = 180
MAX_LAT = 90
ARC = 6371393.0


class Point(object):

    def __init__(self, longitude, latitude):
        """经纬度表示的点"""
        self.__longitude = longitude
        self.__latitude = latitude

    def __repr__(self):
        return 'Point({}, {})'.format(self.longitude, self.latitude)

    def is_in_earth(self):
        return (-MAX_LNG < self.__longitude < MAX_LNG) and (-MAX_LAT < self.__latitude < MAX_LAT)

    @property
    def longitude(self):
        return self.__longitude

    @property
    def latitude(self):
        return self.__latitude


class Circular(object):

    def __init__(self, point, radius):
        """圆形围栏"""
        self.__point = point
        self.__radius = radius
        self.__circumscribed_rectangle = None

    def __repr__(self):
        return 'Circular({}, {})'.format(self.point, self.radius)

    def __contains__(self, point):
        dx = point.longitude - self.__point.longitude
        dy = point.latitude - self.__point.latitude
        b = (point.latitude + self.__point.latitude) / 2.0
        lx = radians(dx) * ARC * cos(radians(b))
        ly = radians(dy) * ARC
        dist = sqrt(lx * lx + ly * ly)
        return dist <= self.__radius

    @property
    def point(self):
        return self.__point

    @property
    def radius(self):
        return self.__radius

    @property
    def circumscribed_rectangle(self):
        if self.__circumscribed_rectangle is None:
            lat_offset = self.__radius * (360 / (ARC * 2 * pi))
            lng_offset = self.__radius * (360 / (ARC * cos(radians(self.__point.latitude)) * 2 * pi))
            self.__circumscribed_rectangle = Rectangle(
                Point(self.__point.longitude - lng_offset, self.__point.latitude + lat_offset),
                Point(self.__point.longitude + lng_offset, self.__point.latitude - lat_offset),
            )
        return self.__circumscribed_rectangle


class Rectangle(object):

    def __init__(self, left_top_point, right_bottom_point):
        """矩形围栏"""
        self.__ltp = left_top_point
        self.__rbp = right_bottom_point

    def __repr__(self):
        return 'Rectangle({}, {})'.format(self.__ltp, self.__rbp)

    def __contains__(self, point):
        return self.__poi_in_rectangle(
            point.longitude, point.latitude,
            self.__ltp.longitude, self.__ltp.latitude,
            self.__rbp.longitude, self.__rbp.latitude,
        )

    @staticmethod
    def __poi_in_rectangle(longitude, latitude, min_longitude, max_latitude, max_longitude, min_latitude):
        if min_latitude <= latitude <= max_latitude:
            if min_longitude * max_longitude > 0:
                if min_longitude <= longitude <= max_longitude:
                    return True
            else:
                if (abs(min_longitude) + abs(max_longitude)) < MAX_LNG:
                    # 跨0度经线
                    if min_longitude <= longitude <= max_longitude:
                        return True
                else:
                    # 跨+/-180度经线
                    left = max(min_longitude, max_longitude)
                    right = min(min_longitude, max_longitude)
                    if left <= longitude <= MAX_LNG or right <= longitude <= -MAX_LNG:
                        return True
        return False

    @property
    def circumscribed_rectangle(self):
        return self

    @property
    def left_top_point(self):
        return self.__ltp

    @property
    def right_bottom_point(self):
        return self.__rbp


class Polygon(object):

    def __init__(self, *points):
        """多边形围栏"""
        self.__points = tuple(points)
        self.__circumscribed_rectangle = None

    def __repr__(self):
        return 'Polygon({})'.format(', '.join([repr(p) for p in self.__points]))

    def __contains__(self, point):
        if not point.is_in_earth():
            return False
        count = 0
        for index in range(len(self.__points)-1):
            if self.__is_ray_intersects_segment(point, self.__points[index], self.__points[index+1]):
                count += 1
        return count % 2 != 0

    @staticmethod
    def __is_ray_intersects_segment(poi, s_poi, e_poi):
        if s_poi.latitude == e_poi.latitude:
            return False
        if s_poi.latitude > poi.latitude and e_poi.latitude > poi.latitude:
            return False
        if s_poi.latitude < poi.latitude and e_poi.latitude < poi.latitude:
            return False
        if s_poi.latitude == poi.latitude and e_poi.latitude > poi.latitude:
            return False
        if e_poi.latitude == poi.latitude and s_poi.latitude > poi.latitude:
            return False
        if s_poi.longitude < poi.longitude and e_poi.longitude < poi.longitude:
            return False
        xseg = (
            e_poi.longitude
            + (poi.latitude - e_poi.latitude)
            / (s_poi.latitude - e_poi.latitude)
            * (s_poi.longitude - e_poi.longitude)
        )
        if xseg < poi.longitude:
            return False
        return True

    @property
    def circumscribed_rectangle(self):
        if self.__circumscribed_rectangle is None:
            lngs = [p.longitude for p in self.__points]
            lats = [p.latitude for p in self.__points]
            self.__circumscribed_rectangle = Rectangle(
                Point(min(lngs), max(lats)),
                Point(max(lngs), min(lats))
            )
        return self.__circumscribed_rectangle

    @property
    def points(self):
        return self.__points


class FenceRTree(object):
    """R-tree for geography fence"""

    def __init__(self):
        self.db = {}  # R-tree database

    def build(self, fences):
        self.db = {}
        self.__init_root_circumscribed_rectangle(fences)
        self.__init_fence_circumscribed_rectangle_tree(fences, line_num=16)

    def __init_root_circumscribed_rectangle(self, fences):
        coordinate = None
        for fence in fences.values():
            rectangle = fence.circumscribed_rectangle
            if coordinate is None:
                coordinate = [
                    rectangle.left_top_point.longitude, rectangle.left_top_point.latitude,
                    rectangle.right_bottom_point.longitude, rectangle.right_bottom_point.latitude
                ]
            coordinate[0] = min(coordinate[0], rectangle.left_top_point.longitude)
            coordinate[1] = max(coordinate[1], rectangle.left_top_point.latitude)
            coordinate[2] = max(coordinate[2], rectangle.right_bottom_point.longitude)
            coordinate[3] = min(coordinate[3], rectangle.right_bottom_point.latitude)

        if coordinate:
            self.db["rectangle"] = Rectangle(
                Point(coordinate[0], coordinate[1]), Point(coordinate[2], coordinate[3])
            )

    @staticmethod
    def __init_fence_circumscribed_rectangle(p_rect, s_rect):
        coordinate = [
            min([p_rect.left_top_point.longitude, s_rect.left_top_point.longitude]),
            max([p_rect.left_top_point.latitude, s_rect.left_top_point.latitude]),
            max([p_rect.right_bottom_point.longitude, s_rect.right_bottom_point.longitude]),
            min([p_rect.right_bottom_point.latitude, s_rect.right_bottom_point.latitude]),
        ]
        return Rectangle(Point(*coordinate[:2]), Point(*coordinate[2:]))

    def __init_fence_circumscribed_rectangle_from_areas(self, sk, sl, line_num, _ij, areas):
        _area = {}
        for k in range(sk, sk + 2):
            for l in range(sl, sl + 2):
                _kl = k * line_num + l
                if areas.get(_kl):
                    if _area.get(_ij) is None:
                        _area[_ij] = {
                            "rectangle": areas[_kl]["rectangle"],
                            "areas": [areas[_kl]]
                        }
                    else:
                        _area[_ij]["rectangle"] = self.__init_fence_circumscribed_rectangle(
                            _area[_ij]["rectangle"], areas[_kl]["rectangle"]
                        )
                        _area[_ij]["areas"].append(areas[_kl])
        return _area

    def __init_fence_circumscribed_rectangle_recursion(self, line_num, areas):
        _area = {}
        _line_num = int(line_num / 2)
        if _line_num > 1:
            for i in range(_line_num):
                sk = i * 2
                for j in range(_line_num):
                    _ij = i * _line_num + j
                    sl = j * 2
                    _area_ = self.__init_fence_circumscribed_rectangle_from_areas(sk, sl, line_num, _ij, areas)
                    _area.update(_area_)
            areas = _area
            return self.__init_fence_circumscribed_rectangle_recursion(_line_num, areas)
        else:
            return list(areas.values())

    def __init_fence_circumscribed_rectangle_tree(self, fences, line_num=16):
        root_rectangle = self.db.get("rectangle")
        if not root_rectangle:
            return

        lng_range = (root_rectangle.right_bottom_point.longitude - root_rectangle.left_top_point.longitude) / line_num
        lat_range = (root_rectangle.left_top_point.latitude - root_rectangle.right_bottom_point.latitude) / line_num
        leaf_area = {}
        for fence_id, fence in fences.items():
            rectangle = fence.circumscribed_rectangle
            _lng_no = int((rectangle.left_top_point.longitude - root_rectangle.left_top_point.longitude) / lng_range)
            _lat_no = int((rectangle.right_bottom_point.latitude - root_rectangle.right_bottom_point.latitude) / lat_range)
            k = _lat_no * line_num + _lng_no
            if leaf_area.get(k) is None:
                leaf_area[k] = {
                    "fence_ids": [fence_id],
                    "rectangle": rectangle,
                }
            else:
                leaf_area[k]["fence_ids"].append(fence_id)
                leaf_area[k]["rectangle"] = self.__init_fence_circumscribed_rectangle(
                    leaf_area[k]["rectangle"],
                    rectangle
                )

        self.db["areas"] = self.__init_fence_circumscribed_rectangle_recursion(line_num, leaf_area)

    def __check_point_in_fence_rtree(self, point, rectangle, areas):
        _in_fence_ids = []
        if point in rectangle:
            for item in areas:
                _rectangle = item["rectangle"]
                _areas = item.get("areas", [])
                _fence_ids = item.get("fence_ids", [])
                _in_fence_ids.extend(
                    self.__check_point_in_fence_rtree(
                        point,
                        _rectangle,
                        _areas
                    )
                )
                _in_fence_ids.extend(_fence_ids)
        return _in_fence_ids

    def check_point_in_tree(self, point):
        rectangle = self.db.get("rectangle")
        if not rectangle:
            return {}
        areas = self.db.get("areas", [])
        return self.__check_point_in_fence_rtree(point, rectangle, areas)


class GeoFence(object):
    """geography fence with R-tree"""

    def __init__(self):
        self.fences = {}
        self.rtree = FenceRTree()

    def update(self, data):
        """update fence data

        @data: dict, {"fence_id", <fence_object>, ...}
        """
        self.fences.update(data)

    def build(self):
        """build R tree according fence data"""
        self.rtree.build(self.fences)

    def check_point_in_fences(self, point):
        """check if point in fences

        @point: Point
        return dict, {"fence_id": <fence object>, ...}, point in fences
        """
        fences = {}
        fence_ids = self.rtree.check_point_in_tree(point)
        for fence_id in fence_ids:
            if point in self.fences[fence_id]:
                fences.update({fence_id: self.fences[fence_id]})
        return fences


def test():
    gf = GeoFence()

    # 圆形围栏 测试点：50, 120 --> '1'
    gf.update({
        '0': Circular(Point(113.0, 25.31), 1000),  # [1, 1, 113.0, 25.31, 1000],
        '1': Circular(Point(120.0, 50), 1000),  # [1, 1, 120.0, 50, 1000],
        '2': Circular(Point(117.30794, 31.79322), 100),  # [1, 1, 117.30794, 31.79322, 100],  # 合肥
    })

    # 矩形围栏 测试点：15, 15 --> '3'
    gf.update({
        '3': Rectangle(Point(10, 20), Point(20, 10)),  # [2, 1, 10, 20, 20, 10],
        '4': Rectangle(Point(20, 30), Point(30, 20)),  # [2, 1, 20, 30, 30, 20],
    })

    # 多边形围栏 测试点：31.51031807840651, 117.55787896484374 --> '5'  # 巢湖
    gf.update({
        '5': Polygon(*[
            Point(117.37935113281249, 31.711477410980944),
            Point(117.44801568359374, 31.674085343031127),
            Point(117.49745416015624, 31.57585952710738),
            Point(117.61830376953124, 31.639016599620224),
            Point(117.82704400390624, 31.587558513127775),
            Point(117.61281060546874, 31.425982924793466),
            Point(117.44526910156249, 31.484557065120843),
            Point(117.31068658203124, 31.587558513127775),
            Point(117.30519341796874, 31.688109134141744)
        ])
    })

    gf.build()

    print("======= 圆形围栏 测试点：50, 120 --> '1'")
    print('result: ', gf.check_point_in_fences(Point(120, 50)))
    print("======= 矩形围栏 测试点：15, 15 --> '3'")
    print('result: ', gf.check_point_in_fences(Point(10, 20)))
    print("======= 多边形围栏 测试点：31.51031807840651, 117.55787896484374 --> '5'")
    print('result: ', gf.check_point_in_fences(Point(117.55787896484374, 31.51031807840651)))

    print('=======')
    print('fences: ', gf.fences)
    print('=======')
    print('tree data: ', gf.rtree.db)


if __name__ == '__main__':
    # test()
    earth_rectangle = Rectangle(Point(-180, 90), Point(180, -90))
    p = Point(180, 30)
    print(p in earth_rectangle)