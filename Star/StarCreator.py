import numpy as np
import scipy
import math
import Star.Function as Fn
import Star.Vars as Gv
import warnings
import copy

class StarCreator(object):
    def __init__(self, primary, secondary, orbit):
        pass

        self.primary = Object(mass=primary["mass"], synchronicity_parameter=1.0, potential=primary["potential"],
                              phi_steps=primary["phis"], theta_steps=primary["thetas"])
        self.secondary = Object(mass=secondary["mass"], synchronicity_parameter=1.0, potential=secondary["potential"],
                                phi_steps=secondary["phis"], theta_steps=secondary["thetas"])

        self.orbit = Orbit(orbital_period=orbit["period"], eccentricity=orbit["eccentricity"],
                           inclination=orbit["inclination"], argument_of_periastron=orbit["argument_of_periastron"])
        self.binary = Binary(primary=self.primary, secondary=self.secondary, orbit=self.orbit, system="eb")

    def create_model(self):
        zp = True if self.binary.binary_morph != "over-contact" else False
        ca = np.pi / 2.0 if self.binary.binary_morph != "over-contact" else np.pi / 4.0

        model = self.binary.get_3d_model_optimized(t_object="both", actual_distance=1.0,
                                                   critical_angle=ca,
                                                   phi_steps=self.primary.phi_steps,
                                                   theta_steps=self.primary.theta_steps,
                                                   zero_point=zp, homo=True)
        self.binary.vertices = model

        if not Fn.empty(model):
            pn = Geometry.normal_estimation(binary_object=self.binary, actual_distance=self.binary.periastron_distance,
                                            vertices=np.array(model["primary"]), t_object="primary")
            sn = Geometry.normal_estimation(binary_object=self.binary, actual_distance=self.binary.periastron_distance,
                                            vertices=np.array(model["secondary"]), t_object="secondary")

            if self.binary.binary_morph == "over-contact":
                ns = Geometry.normal_estimation(binary_object=self.binary,
                                                actual_distance=self.binary.periastron_distance,
                                                vertices=model["system"], t_object="primary")
                self.binary.norms = {"primary": pn, "secondary": sn, "system": ns}

            else:
                self.binary.norms = {"primary": pn, "secondary": sn, "system": np.concatenate((pn, sn), axis=0)}

    def create_spots(self, meta=None):
        spots, sizes, centers, alt_centers, boundaries = [], [], [], [], []
        spots_meta = meta
        x_separation = None

        # if wuma system, get separation and make additional test to location of each point (if primary
        # spot doesn't intersect with secondary, if does, then such spot will be skiped completly)
        if self.binary.binary_morph == "over-contact":
            x_separation = self.binary.get_separation(actual_distance=self.binary.periastron_distance)[0]

        for meta in spots_meta:
            lon, lat, diameter = meta["lon"], meta["lat"], meta["diameter"]  # lon -> phi, lat -> theta
            stps_azimuthal, stps_radial = meta["steps_azimuthal"], meta["steps_radial"]

            # angle to move in radius of spot
            mov_angle = float(diameter * 0.5) / float(stps_radial)

            # radial vector (for spot point computation)
            # radial_v = Fn.spheric_to_cartesian(vector=np.array([1.0, lon, lat]))
            radial_v = np.array([1.0, lon, lat])  # unit radial vector to center of current spot

            boundary, spot, solution = [], [], False
            args, use = (self.binary.periastron_distance, radial_v[1], radial_v[2]), False

            solution, use = self.solve(args, meta["t_object"], x_separation=x_separation)
            # radial distance to the center of current spot from begining of coordinate system
            spot_center_r = solution

            # # for False testing
            # if meta["diameter"] == 0.2:
            #     use = False

            if not use:
                # in case of spots, every point should be usefull, otherwise skip current spot computation
                self.set_false(x=[alt_centers, boundaries, spots, centers, sizes])
                continue

            center = Fn.spheric_to_cartesian(vector=np.array([solution, radial_v[1], radial_v[2]])).tolist()
            spot.append(center)

            # adaptive azimuthal steps in plane of spot
            n0 = int(stps_azimuthal)

            # we have to obtain distance between center and firts point in first circle of spot
            args, use = (self.binary.periastron_distance, lon, lat + mov_angle), False
            solution, use = self.solve(args, meta["t_object"], x_separation=x_separation)

            # # for False testing
            # if meta["diameter"] == 0.2:
            #     use = False
            if not use:
                # in case of spots, every point should be usefull, otherwise skip current spot computation
                self.set_false(x=[alt_centers, boundaries, spots, centers, sizes])
                continue

            x0 = np.sqrt(spot_center_r**2 + solution**2 - (2.0 * spot_center_r * solution * np.cos(mov_angle)))

            # if problem will occured in while loop down bellow, we have to break current spot computation
            # while will break and next_step will set to False
            next_step = True
            for i in range(0, stps_radial):
                spheric_v = np.array([1.0, lon, lat + (mov_angle * (float(i) + 1.0))])

                if spheric_v[1] < 0:
                    spheric_v[1] += 2.0 * np.pi

                # adaptive steps in plane of spot
                ni = n0 * (float(i) + 1.0)
                # ni = n0 * ((float(i) + 1.0) * x0) / x0

                rot_angle, phi = 2.0 * np.pi / ni, 0.0

                while phi < (2.0 * np.pi - (rot_angle * 0.5)):
                    args, use = (self.binary.periastron_distance, spheric_v[1], spheric_v[2]), False
                    solution, use = self.solve(args, meta["t_object"], x_separation=x_separation)

                    # # for False testing
                    # if meta["diameter"] == 0.2:
                    #     use = False

                    if not use:
                        # in case of spots, every point should be usefull, otherwise skip current spot computation
                        self.set_false(x=[alt_centers, boundaries, spots, centers, sizes])
                        next_step = False
                        break

                    current_point = np.array([solution, spheric_v[1], spheric_v[2]])
                    xyz = Fn.spheric_to_cartesian(vector=current_point)

                    # spot.append([self.x_flip(t_object=meta["t_object"], x=xyz[0],
                    #             star_separation=self.binary.periastron_distance), xyz[1], xyz[2]])

                    spot.append(xyz.tolist())

                    cartesian_v = Fn.spheric_to_cartesian(vector=spheric_v)
                    # vector radial_v is normalized in fucntion arbitrary_rotation
                    rot_v = Fn.arbitrary_rotate(theta=rot_angle, vector=cartesian_v,
                                                omega=Fn.spheric_to_cartesian(radial_v))
                    spheric_v = Fn.cartesian_to_spheric(vector=rot_v, degrees=False)

                    if spheric_v[1] < 0:
                        spheric_v[1] += 2.0 * np.pi
                    phi += rot_angle

                    if i == int(stps_radial) - 1:
                        boundary.append(xyz.tolist())

            if not next_step:
                continue

            # alternative center of spot (obtained from spot boundary)
            # boundary center of mass
            # return np.array([(tri[0] + tri[1] + tri[2]) / 3.0 for tri in faces])

            b_com = sum(np.array(boundary)) / len(boundary)
            spheric_v = Fn.cartesian_to_spheric(vector=b_com, degrees=False)

            # i am not testing use, beacuse if points of entire spot exist, then this center point also has to
            args, use, solution = (self.binary.periastron_distance, spheric_v[1], spheric_v[2]), False, False
            solution, use = self.solve(args, meta["t_object"], x_separation=x_separation)
            current_point = np.array([solution, spheric_v[1], spheric_v[2]])
            alt_center = Fn.spheric_to_cartesian(vector=current_point).tolist()

            spot.insert(0, alt_center)
            # ----------------------------

            sizes.append(max([np.linalg.norm(np.array(alt_center) - np.array(b)) for b in boundary]))
            if meta["t_object"] == "primary":
                alt_centers.append(alt_center)
                boundaries.append(boundary)
                spots.append(spot)
                centers.append(center)
            elif meta["t_object"] == "secondary":
                alt_centers.append([self.x_flip(t_object="secondary", x=alt_center[0],
                                                star_separation=self.binary.periastron_distance),
                                    -alt_center[1], alt_center[2]])

                centers.append([self.x_flip(t_object="secondary", x=center[0],
                                            star_separation=self.binary.periastron_distance),
                                -center[1], center[2]])

                b = [[self.x_flip(t_object="secondary", x=x[0],
                                  star_separation=self.binary.periastron_distance),
                      -x[1], x[2]] for x in boundary]
                boundaries.append(b)

                s = [[self.x_flip(t_object="secondary", x=x[0],
                                  star_separation=self.binary.periastron_distance),
                      -x[1], x[2]] for x in spot]
                spots.append(s)

            else:
                return False

        spots_d = {"primary": [], "secondary": []}
        for meta, i in list(zip(spots_meta, range(0, len(spots_meta)))):
            if not Fn.empty(spots[i]):
                norms = Geometry.normal_estimation(binary_object=self.binary,
                                                   actual_distance=self.binary.periastron_distance,
                                                   vertices=np.array(spots[i]), t_object=meta["t_object"])

                spots_d[meta["t_object"]].append(
                    {"vertices": spots[i], "center": centers[i], "alt_center": alt_centers[i],
                     "size": sizes[i], "boundary": boundaries[i], "norms": norms})
        self.binary.spots = spots_d

    @classmethod
    def x_flip(cls, t_object=None, x=None, star_separation=None):
        if t_object == "primary":
            return x
        if t_object == "secondary":
            return - (x - star_separation)

    @classmethod
    def set_false(cls, x):
        for t in x:
            t.append(False)

    def solve(self, args, t_object=None, x_separation=None):
        # args = actual_distance, phi, theta
        args, use, solution = args, False, False
        try:
            scipy_solve_point = self.primary.polar_radius / 10.0 if t_object == "primary" else self.secondary.polar_radius / 10.0
            potential_fn = self.binary.primary_potential_fn if t_object == "primary" else self.binary.secondary_potential_fn
            solution, _, ier, _ = scipy.optimize.fsolve(potential_fn, scipy_solve_point, full_output=True, args=args)

            if ier == 1 and not np.isnan(solution[0]):
                solution = solution[0]
                if 1 > solution > 0:
                    use = True
        except:
            use = False

        # test if point is not "behind" separation in case of wuma (if not wuma, then x_separation is set to None)
        if x_separation is not None and use == True:
            # x value
            x = Fn.spheric_to_cartesian([solution, args[1], args[2]])[0]
            x = x if t_object == "primary" else self.x_flip(t_object="secondary", x=x,
                                                            star_separation=args[0])

            # also equal, i dont care abot bullshit spots
            if (t_object == "primary" and x >= x_separation) or (t_object == "secondary" and x <= x_separation):
                use = False

        return solution, use

    def get_vertices(self):
        return self.binary.vertices

    def get_norms(self):
        return self.binary.norms

    def get_spots(self):
        return self.binary.spots

    def get_binary_morph(self):
        return self.binary.binary_morph


class Object(object):
    def __init__(
            self,
            mass=None,  # unit: solar mass //
            synchronicity_parameter=None,  # unit: [-]
            potential=None,  # unit: [-]
            phi_steps=10,
            theta_steps=20
    ):
        self.phi_steps = int(phi_steps)
        self.theta_steps = int(theta_steps)

        self.mass = float(mass)
        self.synchronicity_parameter = float(synchronicity_parameter)
        self.potential = float(potential)

        self.polar_radius = None
        self.critical_potential = None


class Binary(object):
    def __init__(
            self,
            primary=None,
            secondary=None,
            orbit=None,
            system=None,  # eb/te
            planete="roche",  # roche/sphere
    ):
        self.vertices = None
        self.norms = None
        self.spots = None

        self.exception = []
        self.planete = planete  # zadefinovane zatial len preto, aby pycharm nepapuloval
        self.primary = primary
        self.secondary = secondary
        self.system = system

        self.init = True
        self.binary_morph = None
        self.relative_semimajor_axis = None

        self.mass_ratio = secondary.mass / primary.mass
        self.invert_mass_ratio = 1.0 / self.mass_ratio

        try:
            if orbit is None:
                self.init = False
                self.exception.append("In class: Binary, function: __init__(), line: " +
                                      str(Fn.lineno()) + ". Missing orbit class on initialisation.")
                raise Exception()

            self.orbit = orbit

            # try:
            self.primary.critical_potential = self.critical_potential(t_object="primary",
                                                                      actual_distance=orbit.get_periastron_distance())
            self.secondary.critical_potential = self.critical_potential(t_object="secondary",
                                                                        actual_distance=orbit.get_periastron_distance())

            if self.primary.synchronicity_parameter == 1.0 and self.secondary.synchronicity_parameter == 1.0 \
                    and orbit.get_eccentricity() == 0.0:
                self.primary.lagrangian_points = self.get_lagrangian_points(
                    actual_distance=1.0,
                    synchronicity_parameter=self.primary.synchronicity_parameter,
                    t_object="primary")
                # toto nastavenie sekundarnych lagrangeovych bodov je tu len pre uplnost, aby tam neostala hodnota None
                self.secondary.lagrangian_points = self.primary.lagrangian_points

                argsi, argsii = (1.0, self.primary.lagrangian_points[1], 0., np.pi / 2.), ()
                # podla pomeru hmotnosti sa za L2 berie ten, ktory je za lahsou zlozkou, cize ak q <= 1, tak sa berie
                # lagrangian_points[2] a inak sa berie lagrangian_points[0]
                if 1 >= self.mass_ratio > 0:
                    argsii = (1.0, self.primary.lagrangian_points[2], 0., np.pi / 2.)
                elif self.mass_ratio > 1:
                    argsii = (1.0, abs(self.primary.lagrangian_points[0]), np.pi, np.pi / 2.)

                potential_inner = abs(self.potential_value(*argsi))
                potential_outer = abs(self.potential_value(*argsii))

                df_potential = potential_inner - potential_outer

                # premenna nastavena len pre self testing, nema realne vyuzitie ako self.
                self.potential_inner = potential_inner
                self.potential_outer = potential_outer
                self.df_potential = df_potential

                # nastavenie premennych pre filling faktory primarnej a sekundarnej zlozky
                self.primary.filling_factor = (potential_inner - self.primary.potential) / df_potential
                self.secondary.filling_factor = (potential_inner - self.secondary.potential) / df_potential

                if ((1 > self.secondary.filling_factor > 0) or (1 > self.primary.filling_factor > 0)) and (
                            self.primary.filling_factor != self.secondary.filling_factor):
                    self.init = False
                    self.exception.append("ValueError: In class: Binary, function: __init__(), line: " + str(
                        Fn.lineno()) + ". Detected over-contact system, but potentials of components don't match.")
                    raise Exception()

                if self.primary.filling_factor > 1 or self.secondary.filling_factor > 1:
                    self.exception.append(
                        "ValueError In class: Binary, function: __init__(), line: " + str(Fn.lineno()) +
                        ". Non-Physical system: primary.filling_factor > 1 or secondary.filling_factor > 1")
                    self.init = False
                    raise Exception()

                if self.primary.filling_factor < 0 and self.secondary.filling_factor < 0:
                    self.binary_morph = "detached"
                elif (self.primary.filling_factor == 0 and self.secondary.filling_factor < 0) or (
                                self.primary.filling_factor < 0 and self.secondary.filling_factor == 0):
                    self.binary_morph = "semi-contact"
                elif 1 > self.primary.filling_factor > 0:
                    self.binary_morph = "over-contact"
                elif self.primary.filling_factor > 1 or self.secondary.filling_factor > 1:
                    self.init = False
                    self.exception.append("Info: Open system. Nothing more to do.")
            else:
                self.binary_morph = "WARNING, NOT CIRCULAR ORBIT"

            self.periastron_distance = orbit.get_periastron_distance()
            self.primary.polar_radius = self.get_polar_radius(t_object="primary",
                                                              actual_distance=self.periastron_distance)
            self.secondary.polar_radius = self.get_polar_radius(t_object="secondary",
                                                                actual_distance=self.periastron_distance)

            if not self.primary.polar_radius or not self.secondary.polar_radius:
                self.exception.append("")
                raise Exception()

            self.primary.backward_radius = self.get_backward_radius(t_object="primary",
                                                                    actual_distance=self.periastron_distance)
            self.secondary.backward_radius = self.get_backward_radius(t_object="secondary",
                                                                      actual_distance=self.periastron_distance)

            # relativna dlzka hlavnej polosi sa nastavi pre instanciu Binary aj pre instanciu Orbit
            self.relative_semimajor_axis = (((((self.orbit.orbital_period * 86400.0) ** 2) * (
                Gv.G_CONSTANT * Gv.SOLAR_MASS * (self.primary.mass + self.secondary.mass))) / (
                                                 4.0 * np.pi ** 2)) ** (
                                                1.0 / 3.0)) / Gv.SOLAR_RADIUS  # in Gv.SOLAR_RADIUS unit
            self.orbit.relative_semimajor_axis = self.relative_semimajor_axis

        except:
            self.init = False
            self.exception.append("Error ocured in initialisation of binary system.")

    def critical_potential(self, t_object="primary", actual_distance=None):
        solution = None

        if t_object == "primary":
            args = (actual_distance, self.primary.synchronicity_parameter)
            solution = scipy.optimize.newton(self.primary_potential_derivation_x, 0.001, args=args)
        if t_object == "secondary":
            args = (actual_distance, self.secondary.synchronicity_parameter)
            solution = scipy.optimize.newton(self.secondary_potential_derivation_x, 0.001, args=args)
        if not np.isnan(solution):
            if t_object == "primary":
                args = (actual_distance, solution, 0.0, np.pi / 2.)
                return abs(self.potential_value(*args))
            if t_object == "secondary":
                args = (actual_distance, actual_distance - solution, 0.0, np.pi / 2.)
                return abs(self.inverted_potential_value(*args))
        else:
            self.exception.append("ValueError: In class: Binary, function: critical_potential(), line: " + str(
                Fn.lineno()) + ". Wrong value has been encoutered.")
            return False

    def wc_potential_fn(self, radius_p, *args):
        actual_distance, phi, x = args
        F = 1.0
        potential_omega = self.primary.potential

        block_a = (1.0 / np.sqrt(radius_p**2 + x**2))
        block_b = self.mass_ratio * ((1.0 / (np.sqrt((x - actual_distance)**2 + radius_p**2))) -
                                     (x / actual_distance**2))
        block_c = 0.5 * (1.0 + self.mass_ratio) * (F**2) * (x**2 + (radius_p**2 * np.cos(phi)**2))
        return block_a + block_b + block_c - potential_omega

    def primary_potential_fn(self, radius, *args):
        # variables
        actual_distance, phi, theta = args
        # /variables
        # block of function
        block_a = (1 / radius)
        block_b = (self.mass_ratio / (np.sqrt(np.power(actual_distance, 2) + np.power(radius, 2) - (
            2. * radius * (np.cos(phi) * np.sin(theta)) * actual_distance))))
        block_c = ((self.mass_ratio * radius * (np.cos(phi) * np.sin(theta))) / (np.power(actual_distance, 2)))
        block_d = (
            0.5 * np.power(self.primary.synchronicity_parameter, 2) * (1 + self.mass_ratio) * np.power(radius, 2) * (
                1 - np.power(np.cos(theta), 2)))
        # /block of function
        return (block_a + block_b - block_c + block_d) - self.primary.potential

    def secondary_potential_fn(self, radius, *args):
        # variables
        actual_distance, phi, theta = args
        # /variables
        # block of function
        block_a = (1. / radius)
        block_b = (self.invert_mass_ratio / (np.sqrt(np.power(actual_distance, 2) + np.power(radius, 2) - (
            2 * radius * (np.cos(phi) * np.sin(theta)) * actual_distance))))
        block_c = ((self.invert_mass_ratio * radius * (np.cos(phi) * np.sin(theta))) / (np.power(actual_distance, 2)))
        block_d = (
            0.5 * np.power(self.secondary.synchronicity_parameter, 2) * (1 + self.invert_mass_ratio) * np.power(
                radius,
                2) * (
                1 - np.power(np.cos(theta), 2)))
        # /block of function
        inverse_potential = (block_a + block_b - block_c + block_d) / self.invert_mass_ratio + (
            0.5 * ((self.invert_mass_ratio - 1) / self.invert_mass_ratio))
        return inverse_potential - self.secondary.potential

    def polar_potential_fn(self, radius, *args):
        # variables
        actual_distance, t_object = args
        # /variabreturn - (block_a + block_b - block_c + block_d)les
        if t_object == "primary":
            return (1. / radius) + (self.mass_ratio * (
                (1 / (np.sqrt(np.power(actual_distance, 2) + np.power(radius, 2)))))) - self.primary.potential
        if t_object == "secondary":
            return (((1. / radius) + (
                self.invert_mass_ratio * (1 / (np.sqrt(np.power(actual_distance, 2) + np.power(radius, 2)))))) / (
                        self.invert_mass_ratio) + (
                        0.5 * ((self.invert_mass_ratio - 1) / self.invert_mass_ratio))) - self.secondary.potential

    def get_polar_radius(self, t_object="primary", actual_distance=None):

        # premenna scipy_solve_point sluzi ako startovaci bod pre numercky vypocet implicitnej rovnice pre polarny
        # polomer prvykrat pri vyvoji bola nastavena na nieco okolo 0.01 no pri strasne malych zlozkach, ktorych
        # polomer je pdo touto hodnotou numerika neskonvergovala, tak musel byt nastavena na hodnotu vyrazne mensiu,
        # vsetky dalsie vypocty po polarnom polomere (uz ked je znamy) budu musiet byt nastavene na nejaku 1/10 hodnoty
        # polomeru pre danu zlozku systemu
        ret, scipy_solve_point, solution = False, 1e-20, "NaN"
        np.seterr(divide='raise', invalid='raise')
        try:
            args = (actual_distance, t_object)
            solution, info, ier, msg = scipy.optimize.fsolve(self.polar_potential_fn, scipy_solve_point,
                                                             full_output=True, args=args)
            if ier == 1 and not np.isnan(solution[0]):
                solution = solution[0]
                # tu sa kontroluje, ci je spocitany polomer vacsi ako 0 a mesni ako 1
                # urcite nemoze byt mensi ako 0, a nepredpoklada sa, ze presiahen hodnotu separacie 1 pri zlozkach,
                # takyto system asi nemoze nijak vzniknut;
                # ak je hodnota ina, pojde urcite o chybu numeriky
                if 1 > solution > 0:
                    ret = True
        except:
            ret = False
        np.seterr(divide='print', invalid='print')
        if ret:
            return solution
        return False

    def get_backward_radius(self, t_object="primary", actual_distance=None):
        args, ret, ier, solution = (actual_distance, np.pi, np.pi / 2.0), False, 10, 0.0
        try:
            if t_object == "primary":
                scipy_solve_point = self.primary.polar_radius / 10.0
                solution, info, ier, msg = scipy.optimize.fsolve(self.primary_potential_fn,
                                                                 scipy_solve_point,
                                                                 full_output=True,
                                                                 args=args)

            elif t_object == "secondary":
                scipy_solve_point = self.secondary.polar_radius / 10.0
                solution, info, ier, msg = scipy.optimize.fsolve(self.secondary_potential_fn,
                                                                 scipy_solve_point,
                                                                 full_output=True, args=args, )

            if ier == 1 and not np.isnan(solution[0]):
                solution = solution[0]
                # tu sa kontroluje, ci je spocitany polomer vacsi ako 0 a mesni ako 1
                # urcite nemoze byt mensi ako 0, a nepredpoklada sa, ze presiahen hodnotu separacie 1 pri zlozkach,
                # takyto system ais nemoze nijak vzniknut; ak je hodnota ina, pojde urcite o chybu numeriky
                if 1 > solution > 0:
                    ret = True
        except:
            ret = False
        if ret:
            return solution
        return False

    def primary_potential_derivation_x(self, x, *args):
        actual_distance, synchronicity_parameter = args
        r_sqr, rw_sqr = x ** 2, (actual_distance - x) ** 2
        return - (x / r_sqr ** (3. / 2.)) + (
            (self.mass_ratio * (actual_distance - x)) / rw_sqr ** (3. / 2.)) + synchronicity_parameter ** 2 * (
            self.mass_ratio + 1) * x - self.mass_ratio / actual_distance ** 2

    def secondary_potential_derivation_x(self, x, *args):
        actual_distance, synchronicity_parameter = args
        r_sqr, rw_sqr = x ** 2, (actual_distance - x) ** 2
        return - (x / r_sqr ** (3. / 2.)) + (
            (self.mass_ratio * (actual_distance - x)) / rw_sqr ** (3. / 2.)) - synchronicity_parameter ** 2 * (
            self.mass_ratio + 1) * (1 - x) + (1. / actual_distance ** 2)

    def potential_value(self, *args):
        # variables
        actual_distance, radius, phi, theta = args

        # block of functions
        block_a = (1. / radius)
        block_b = (self.mass_ratio / (np.sqrt(np.power(actual_distance, 2) + np.power(radius, 2) - (
            2. * radius * (np.cos(phi) * np.sin(theta)) * actual_distance))))
        block_c = ((self.mass_ratio * radius * (np.cos(phi) * np.sin(theta))) / (np.power(actual_distance, 2)))
        block_d = (
            0.5 * np.power(self.primary.synchronicity_parameter, 2) * (1 + self.mass_ratio) * np.power(radius, 2) * (
                1 - np.power(np.cos(theta), 2)))
        return - (block_a + block_b - block_c + block_d)

    def inverted_potential_value(self, *args):
        # variables
        actual_distance, radius, phi, theta = args
        # block of functions
        block_a = (1. / radius)
        block_b = (self.invert_mass_ratio / (np.sqrt(np.power(actual_distance, 2) + np.power(radius, 2) - (
            2 * radius * (np.cos(phi) * np.sin(theta)) * actual_distance))))
        block_c = ((self.invert_mass_ratio * radius * (np.cos(phi) * np.sin(theta))) / (np.power(actual_distance, 2)))
        block_d = (
            0.5 * np.power(self.secondary.synchronicity_parameter, 2) * (1 + self.invert_mass_ratio) * np.power(
                radius,
                2) * (
                1 - np.power(np.cos(theta), 2)))
        inverse_potential = (block_a + block_b - block_c + block_d) / self.invert_mass_ratio + (
            0.5 * ((self.invert_mass_ratio - 1) / self.invert_mass_ratio))
        return - inverse_potential

    def get_lagrangian_points(self, solve_step=0.1, actual_distance=None, synchronicity_parameter=1.0,
                              t_object='primary'):

        # ta funkcia funguje tak, ze zacne ratat numericky derivaciu potencialu po x - ovej osi
        # problem je v tom, ze to z nejakej polohy zrata L2, z inej L1 a z inej L3 v zmysle ineho startovacieho bodu
        # z toho dovodu musim prebehnut po osi s nejakym krokom a nechat vyriesit numeriku v kazdom bode
        # dalsi problem je v tom, ze kedze sa jedna o numeriku, tie doratane hodnoty nie su uplne presne
        # teda sa nedaju jednoducho porovnat stylom L_(n) == L_(n-1), pretoze by sa nemuseli rovnat aj keby to boli
        # tie iste lagrangeove body, preto to treba nejak rozumen zaokruhil a tak to porovnavat a "presne" hodnoty
        # ukladat inde
        args, end = (actual_distance, synchronicity_parameter), int(((actual_distance * 6.) / solve_step))
        points, lagrange = [], []
        scipy_solve_point, decimal = - (actual_distance * 3.) + solve_step, int(len(str(solve_step).split('.')[1]))
        derivation_value, ier = None, np.inf

        for i in range(0, end):
            # toto sa vykonava aby sa zistilo, ci nedochadza k deleniu nulou pri danej polohe riesenia rovnice, ak ano,
            # tak sa to proste preskoci
            try:
                np.seterr(divide='raise', invalid='raise')
                if t_object == 'primary':
                    self.primary_potential_derivation_x(round(scipy_solve_point, decimal), *args)
                if t_object == 'secondary':
                    self.secondary_potential_derivation_x(round(scipy_solve_point, decimal), *args)
                np.seterr(divide='print', invalid='print')
                pass
            except:
                scipy_solve_point += solve_step
                continue

            try:
                if t_object == 'primary':
                    solution, info, ier, msg = scipy.optimize.fsolve(self.primary_potential_derivation_x,
                                                                     scipy_solve_point, full_output=True, args=args)
                if t_object == 'secondary':
                    solution, info, ier, msg = scipy.optimize.fsolve(self.secondary_potential_derivation_x,
                                                                     scipy_solve_point, full_output=True, args=args)

                if ier == 1:
                    if round(solution[0], 5) not in points:
                        try:
                            if t_object == 'primary': derivation_value = abs(
                                round(self.primary_potential_derivation_x(solution[0], *args), 4))
                            if t_object == 'secondary': derivation_value = abs(
                                round(self.secondary_potential_derivation_x(solution[0], *args), 4))
                            if derivation_value == 0:
                                use = True
                            else:
                                use = False
                        except:
                            use = False
                        if use:
                            points.append(round(solution[0], 5))
                            lagrange.append(solution[0])
            except:
                scipy_solve_point += solve_step
                continue

            scipy_solve_point += solve_step
        return np.array(sorted(lagrange))

    def get_exception(self):
        return self.exception

    def get_3d_model_optimized(self, phi_steps=10, theta_steps=10, t_object="both", actual_distance=None,
                               critical_angle=np.pi / 4., additional_angle=np.pi / 180., zero_point=True, homo=True):

        # inicializacie
        primary_phi_steps, primary_theta_steps = None, None
        secondary_phi_steps, secondary_theta_steps = None, None
        use = False
        z_point, rotation_angle, transform = None, None, None

        # first adaptive angle (pre vypocet separacnej hranice)
        fat = {"primary": None, "secondary": None}

        # phi_steps, pocet azimutalnych krokov, jedna sa o mierne zavadzajuce pomenovanie, pretoze sa jedna o
        # pocet krokov
        # o kolko sa bude rotovat po kruznici okolo x-ovej osi
        phi_steps = phi_steps
        model, partial, total = {'primary': 'empty', 'secondary': 'empty',
                                 'system': 'empty', "separation": "empty"}, [], [],
        point_coordinate = np.arange(3, dtype=np.float)
        ppr, spr = self.primary.polar_radius, self.secondary.polar_radius

        if actual_distance is None:
            print("ValueError: In class: Binary, function: get_3d_model_testing(), line: " + str(
                Fn.lineno()) + ". Variable `actual_distance` is set to None.")
            return False

        # vypocet bodov na povrchu hviezdy prebieha z polohy phi = np.pi, theta = np.pi/2., to znamena z jej zadnej
        # casti
        # vypocet po zmene polarneho kroku po kruznici okolo x-ovej osi, ktora sa ziska transformaciou
        # critical_angle sluzi pre zhustenie vypoctu na hrdle wuma systemov, ked polarny uhol narazi na hodnotu
        # np.pi + critical_angle
        # (je to np.pi + critical_angle pretoze ak vypocet zacina na np.pi/2 a pridava sa polarny uhol, tak np.pi
        # dosiahne na hranici
        # tretieho a stvrteho kvadrantu z rovine z-x a dalej teda presahuje np.pi, co nie je problem, pretoze cos
        # je parna funkcia
        # a teada nadobuda rovnaku hodnotu ako pre np.pi + nejaky_uhol ako pre np.pi - nejaky uhol)
        # tu sa skontroluje, ci kriticky uhol nie je vacsi ako pi / 2
        # kriticky uhol sa pocita od 180 stupnov polarneho uhla, vid obrazok
        # oA               |\              Bo
        #  o               | \             o
        #    o             |  \          o
        #      o           | K \       o
        #          o       |    \  o
        #               o  |  o
        # K oznacuje velkost kritickeho uhla
        # ak je K = np.pi/2.0, tak sa pocita klasicky polarny cyklus z polohy A az do polohy B
        if critical_angle > np.pi / 2.: critical_angle = np.pi / 2.

        if t_object == "both":
            j_range_bottom, j_range_top = 0, 2

            if homo:
                # homogenny model, snahavytvorit model tak, aby sekundarna aj primarna zlozka mali homogenne
                # plosne elementy
                # ak je premenna homo nastavena na True, tak sa pocet krokov v polarnom a azimutalnom smere
                # nastavi nasledovne
                # 1) ak je polarny polomer primarnej zlozky mensi ako sekundarny, tak na tejto mensje zlozke
                # (primarnej) bude
                #    zadany pocet polarnych a azimutalnych krokov zo vstupu
                #    na sekundarnej zlozke sa sa zvysi pocet azimutalnych a polarnych krokov faktorom spocitanym
                # z pomeru ich vzajomnych
                #    polomerov
                # 2) ak je to naopak, tak sa prepocita pocet krokov na primarnej zlozke
                # myslienka je udrzat homogenne pokrytie a to tak, aby bolo co najviac bodov na hviezdach
                if ppr <= spr:
                    primary_phi_steps, primary_theta_steps = phi_steps, theta_steps
                    # zaokruhluje sa pocet krokov hore, radsej hustejsie pokrytie ako redsie
                    secondary_phi_steps, secondary_theta_steps = math.ceil((spr / ppr) * phi_steps), math.ceil(
                        (spr / ppr) * theta_steps)
                elif ppr > spr:
                    secondary_phi_steps, secondary_theta_steps = phi_steps, theta_steps
                    primary_phi_steps, primary_theta_steps = math.ceil((ppr / spr) * phi_steps), math.ceil(
                        (ppr / spr) * theta_steps)
            else:
                # ak je premenna homo nastavena na False, tak sa neudrziava homogenita na zlozkach, ale na obe zlozky
                # sa nastavi pocet azimutalnych a polarnych krokov zo vstupu, prakticky teda mensia zlozka bude pokryta
                # hustejsie ako zlozka, ktora ma vasi polomer
                primary_phi_steps, primary_theta_steps = phi_steps, theta_steps
                secondary_phi_steps, secondary_theta_steps = phi_steps, theta_steps
        # v pripade, ze je snaha spocitat len jednu zo zloziek, tak sa nastavia parametre na danu zlozku
        # rovnako tak rozsah vonkajsieho for cyklu sa nastavi na jedno zbehnutie zodpovedajuce danej zlozke
        elif t_object == "primary":
            j_range_bottom, j_range_top = 0, 1
            primary_phi_steps, primary_theta_steps = phi_steps, theta_steps
        elif t_object == "secondary":
            j_range_bottom, j_range_top = 1, 2
            secondary_phi_steps, secondary_theta_steps = phi_steps, theta_steps
        else:
            print("ValueError: In class: Binary, function: get_3d_model_testing(), line: " + str(
                Fn.lineno()) + ". Incorrect `t_object` parameter value.")
            return False

        # for cyklus prebehne cez objekty (primarnu a sekundarnu zlozku)
        for obj in range(j_range_bottom, j_range_top):
            # current_phi_steps a current_theta_steps je pocet azimutalnych a polarnych krokov pre danu hodnotu cyklus
            # teda bude pre primarnu alebo pre sekundarnu zlozku
            if obj == 0:
                current_phi_steps, current_theta_steps = primary_phi_steps, primary_theta_steps
            elif obj == 1:
                current_phi_steps, current_theta_steps = secondary_phi_steps, secondary_theta_steps

            # startovaci uhol theta
            # dlzka poalrneho kroku
            # partial je pole, do ktoreho sa budu ukladat kooordinaty bodov
            theta, theta_step_length, partial, equatorial, meridional = np.pi / 2., np.pi / current_theta_steps, [], [], []

            # posun v polarnom uhle cez rovinu x-z
            #   zo
            #   | o
            #   |  o
            #   |    o
            #   |      o
            #   |         o ->
            #   _________________x
            #
            # postupne tak ako je to zobrazene na obrazku

            for rot in range(0, int(current_theta_steps)):
                # vytvotenie prazdneho vektora pre sfericke a pre kartezianske suradnice
                vector_spheric, vector_xyz = np.arange(3, dtype=np.float), np.arange(3, dtype=np.float)
                # naplnenie spferickeho vektora vychdzou hodnotou pre vypocet na dalsej a dalsej kruznici
                # (dana novou theta)
                # toto osetrenie cez if je tu len pre uplnost, pretoze povodne to bolo takto, ze tam bolo len
                # vector_spheric[1], vector_spheric[2] = 1., np.pi, theta
                # tu je problem, pretoze je divne mat ph = 180 a s napr. theta = 190, su to trochu protichodne uhly
                # funkcia pre transformacie sa s tym vysporiadala nejak tak, ze to fungovalo, ale nepacila sa mi
                # korektnost
                # ked je raz vo sferickych, theta e <0, np.pi>, phi e <0, 2* np.pi>, tak nech je to dodrzane
                if theta <= np.pi:
                    vector_spheric[0], vector_spheric[1], vector_spheric[2] = 1., np.pi, theta
                elif theta > np.pi:
                    vector_spheric[0], vector_spheric[1], vector_spheric[2] = 1.0, 0.0, (2.0 * np.pi) - theta

                # nastavenie hlavneho polomeru pre dalsie vypocty podla momentalneho objektu
                if obj == 0:
                    star_radius = ppr
                elif obj == 1:
                    star_radius = spr

                # cast, ked je polarny uhol mensi ako np.pi + critical_angle
                # polarny uhol bezi stale ako np.pi == 90, 90 + current_theta_steps, 90 + 2 * current_theta_steps
                # ... atd
                # keby nebol kriticky uhol, tak dosiahne az po theta = 270 stupnov, cize prejde celu akoby spodnu
                # vetvu roviny x-z
                # ak je nastaveny kriticky uhol na nejaku hodnotu, tak klasicke ekvidistancne krokovanie pojde len
                # pokial plati
                # podmienka theta <= np.pi + critical_angle, potom sa spravia iste kontroly (ak tato podmienka ine
                # je splnena)
                if theta <= np.pi + critical_angle:
                    # vypocet polomeru v theta
                    # z
                    # |o     |      |
                    # |  o   |      |
                    # |    o | r    | star_radius
                    # |       o     |
                    # |          o  |
                    # |____________ o _________ x
                    # neviem preco tu bolo povodne zaokruhlenie
                    # r = round( star_radius * np.sin( theta - np.pi / 2. ), 4 )

                    r = star_radius * np.sin(theta - np.pi / 2.0)
                    # pocet krokov na kruznici v r (resp v theta), prepocitane pomerovo k polomeru danej mensej kruznice
                    # je to opodmeinkovane tak, ze ak sa jedna o prvy bod, ked sa zacina riesit, tak urcite treba jeden
                    # bod
                    # ale kedze je tam polomer r = 0, tak by aj pocet bodov bol nula, preto tam treba nastavit na
                    # tvrdo "1"
                    phi_points = 1 if theta == np.pi / 2.0 else math.ceil((r / star_radius) * current_phi_steps)
                    # pocet transformacnych krokov po kruznici
                    # opat opodmienakovane, ak sa jedna o prvy bod, tak je pocet transformacnych korkov 1 (keby bola 0,
                    # tak sa cyklus nevykona ani raz)
                    # v ostatnych pripadoch polarneho uhla theta treba + 1 transformacny krok
                    transform_steps = int(phi_points) if theta == np.pi / 2. else int(phi_points) + 1
                    # dlzka kroku, ktory sa bude pouzivat pri transformacii po kruznici
                    # musi sa stale menit, lebo sa meni polomer kruznice a stale sa meni phi_points
                    # nepodstatny je v principe prva dlzka, pretoze je to v zadnej casti hviezdy, phi = np.pi,
                    # theta = np.pi / 2.0
                    phi_step_length = (np.pi / 2.) / phi_points
                    # rotacny uhol je to iste ako phi_step_length, je to tu len znova priradene, kvoli nazvu premennej,
                    # ak by som to nechcel
                    # pre ilustracne ucely a lepsi prehlad v kode, treba to spravit tak, ze "rotation_angle =
                    # ( np.pi / 2. ) / phi_points"
                    # a riadok "rotation_angle = phi_step_length" odstranit
                    rotation_angle = phi_step_length

                ##############################################

                #   TU PRIDE TA KOMPLIKOVANA ADAPTIVNA CAST

                ##############################################

                else:
                    # podmineka theta <= np.pi + critical_angle: nie je splnena ale na druhej strane je theta
                    # mensia ako 270 stupnov,
                    # teda pokial sa theta uhol neprejde cely np.pi od np.pi / 2.0
                    #                           O                                      |z
                    # - np.pi/2. o               x                o  (3./2.) * np.pi   |
                    #              o                           o                      O|_____x
                    #                 o                     o
                    #                     o            o
                    #                          o   o
                    if theta < (3. / 2.) * np.pi:
                        # testovaci polarny uhol zacina z pozicie np.pi / 2.0 (takze z tej, v ktorej moze nastat
                        # problem pri wuma systemoch, na krku
                        # respektive v krku, kde ziaden bod povrchu nie je)
                        testing_theta = np.pi / 2.0
                        # podmienka while cyklu zebezpeci beh len po kriticky uhol, aby to nebezalo do nekonecna
                        while True and testing_theta < (np.pi - critical_angle):
                            # argumenty pre vypocet bodu na povrchu hviezdy
                            args = (actual_distance, 0.0, float(testing_theta))
                            # vypocet bodu na povrchu podobne ako v kode vyssie
                            # pre rozsiahlejsi a lepsi vystup pouzita funkcia fsolve
                            try:
                                if obj == 0:
                                    scipy_solve_point = self.primary.polar_radius / 10.0
                                    solution, info, ier, msg = scipy.optimize.fsolve(self.primary_potential_fn,
                                                                                     scipy_solve_point,
                                                                                     full_output=True,
                                                                                     args=args)
                                elif obj == 1:
                                    scipy_solve_point = self.secondary.polar_radius / 10.0
                                    solution, info, ier, msg = scipy.optimize.fsolve(self.secondary_potential_fn,
                                                                                     scipy_solve_point,
                                                                                     full_output=True, args=args, )
                                if ier == 1 and not np.isnan(solution[0]):
                                    # ak sa naslo korektne riesenie, tak sa ukonci tento cyklus
                                    if 1 > solution[0] > 0:
                                        break
                                # posun v testovacom uhle o 1 stupen (mozno by bolo lepsie nastavit mensi krok,
                                # ale to by zas spomalilo algoritmus)
                                # ak sa nenaslo korektne riesenie, tak cyklus pokracuje
                                testing_theta += np.pi / 180.0
                            except:
                                # ak sa nenaslo korektne riesenie, tak cyklus pokracuje
                                testing_theta += np.pi / 180.0
                                continue

                        # prvy uhol theta pre adaptivny vypocet sa nastavi ten, kde sa prvykrat v cykle vyssie
                        # podarilo najst korektne riesenie
                        # samotny adaptivny uhol zacne tam, kde sa skoncilo ratanie klasickym cyklom na kritickom uhle
                        first_adaptive_theta, theta_adaptive, stop = testing_theta, ((2. * np.pi) - theta), False  #
                        # ak sa zisti, ze uz na theta = np.pi / 2.0 ja koretkne riesenie, tak sa zmeni kriticky
                        # uhol na np.pi / 2.0, to znamena,
                        # ze kalsicky cyklus sa zastavi az na 270 stupnoch
                        if first_adaptive_theta == np.pi / 2.0: critical_angle = np.pi / 2.0  # change critical_
                        # angle if type is EA or EB

                        # cyklus bezi do nastavenia premennej stop na True
                        # este je to opodmienkovane, aby nezbehol, ak sa vyssie nastavi kriticky uhol na np.pi / 2.0

                        # CYKLUS WHILE POLARNEHO UHLA ADAPTIVNEJ CASTI
                        while not stop and first_adaptive_theta != np.pi / 2.0:
                            # spocitanie polomeru v danom uhle theta (theta_adaptive)
                            r = abs(star_radius * np.sin(theta_adaptive - (np.pi / 2.)))

                            # nastavenie premennych pre zrotovanie po kruznici a vypocet bodov na povrchu hviezdy
                            # na tejto kruznici
                            phi_points = math.ceil((r / star_radius) * current_phi_steps)
                            transform_steps = int(phi_points) + 1
                            phi_step_length = (np.pi / 2.) / phi_points
                            rotation_angle = phi_step_length

                            # ked v cykle uz vybehne adaptivny uhol mimo (do oblasti, kde sa neda pocitat (hrdlo WUMa),
                            # tak sa este od kritickeho (prveho adaptivneho) uhla spocita jedna hodnota mierne
                            # posunuta o maly uhol, vlastne hranica akoby)
                            # je to na to, aby sa co najviac priblizilo vypoctom k hranici hrdla wuma systemov
                            # premenan stop sa nastavi na True, teda po dokonceni posledneho cyklu po kruznici sa
                            # uz viac nevykona
                            if theta_adaptive <= first_adaptive_theta:
                                theta_adaptive, stop = first_adaptive_theta + additional_angle, True

                            vector_spheric, vector_xyz = np.arange(3, dtype=np.float), np.arange(3, dtype=np.float)
                            # tu je na rozdiel od klasickeho cyklu poalrneho uhla spravny polarny uhol, nebezi cez
                            # 180 stupnov
                            vector_spheric[0], vector_spheric[1], vector_spheric[2] = 1., 0.0, theta_adaptive

                            # klasicky cyklus po kruznici v polohe "r"
                            # jeho popis je nizzsie pri klasickom cykle, ked sa nerata s kritickym uhlom
                            for transform_adaptive in range(0, transform_steps):
                                args, use = (actual_distance, vector_spheric[1], vector_spheric[2]), False
                                try:
                                    # pre istotu je v tejto casti pouzita funkcia pre riesenie korenov fsolve,
                                    # pretoze newton moze zratat korektnu hodnotu,
                                    # ktora, ale nelezi na hviezde ale niekde vo vonakjsej casti hillovho priestoru
                                    if obj == 0:
                                        scipy_solve_point = self.primary.polar_radius / 10.
                                        solution, info, ier, msg = scipy.optimize.fsolve(self.primary_potential_fn,
                                                                                         scipy_solve_point,
                                                                                         full_output=True, args=args)
                                    elif obj == 1:
                                        scipy_solve_point = self.secondary.polar_radius / 10.
                                        solution, info, ier, msg = scipy.optimize.fsolve(self.secondary_potential_fn,
                                                                                         scipy_solve_point,
                                                                                         full_output=True, args=args)
                                    if ier == 1 and not np.isnan(solution[0]):
                                        solution = solution[0]
                                        if actual_distance >= solution >= 0: use = True
                                except:
                                    use = False
                                if use:
                                    point_coordinate[0], point_coordinate[1], point_coordinate[2] = solution, \
                                                                                                    vector_spheric[1], \
                                                                                                    vector_spheric[2]
                                    xyz = Fn.spheric_to_cartesian(point_coordinate)
                                    if obj == 0:
                                        x = xyz[0]
                                    elif obj == 1:
                                        x = - (xyz[0] - actual_distance)
                                    y, z = xyz[1], xyz[2]
                                    if transform_adaptive == transform_steps - 1:
                                        equatorial.append([x, y, z])
                                    elif transform_adaptive == 0:
                                        meridional.append([x, y, z])
                                    else:
                                        partial.append([x, y, z])
                                vector_xyz = Fn.spheric_to_cartesian(vector_spheric)
                                rotate = Fn.rotate(vector=vector_xyz, angle=rotation_angle)
                                vector_spheric = Fn.cartesian_to_spheric(rotate)
                                if vector_spheric[1] < 0: vector_spheric[1] += (2. * np.pi)
                            # theta adaptive sa defaultne prvykrat nastavil na  2.0 * np.pi - theta, cize sa
                            # znormalizovala
                            # jeho pozicia na menej ako 180 stupnov, cize ak mal napr. polarny uhol theta hodnotu
                            # 225 stupnov
                            # (samozrjeme vsetko v radianoch), tak sa jeho hodnota zmenila na 135 stupnov, v tomto cykle
                            # polarneho uhla (je nim ten while vyssie) sa teda odpocitava z hodnoty uhla, aby sme sa
                            # pohybovali s uhlom
                            # v smere tak ako v povodnom for cykle, teda smerom od zapornej casti z osi ku kladnej
                            # casti z pohladu
                            # roviny x-z
                            # to preco sa odpocitava taky zvlastny krok??? preto lebo ked som to skusal, tak toto sa
                            # dalo pouzit,
                            # mozno by sa dalo este aj nieco logaritmicke, ale toto funguje
                            theta_adaptive -= np.sin(theta_adaptive) * theta_step_length / 1.5

                        # toto je osetrenie, aby v pripade, ze naozaj bolo pouzite adaptivne pocitanie nedoslo k
                        # dokonceniu klasickeho for cyklus
                        # pre polarny uhol, pretoze k dokonceniu vypoctu povrchu hviezdnej zlozky doslo adaptvne
                        # tato cast je uz mimo adaptivnu, preto po break-u tu, dojde k ukonceniu polarneho for cyklu
                        if first_adaptive_theta > np.pi / 2.:
                            break

                # cyklus po kruznici v polomere "r" (t.j. v uhle theta)
                for transform in range(0, transform_steps):
                    # argumenty pre funkciu potencialu
                    # premenna use pre rozhodnutie pouzitia vysledku numerickeho vypoctu (nemusi stale skonvergovat)
                    args, use = (actual_distance, vector_spheric[1], vector_spheric[2]), False
                    # pokus o vyriesenie implicitnej rovnice (pokus, preto je to v try)
                    try:
                        # pouzije sa prislusny potencial podla zlozky
                        if obj == 0:
                            scipy_solve_point = self.primary.polar_radius / 10.
                            solution = scipy.optimize.newton(self.primary_potential_fn, scipy_solve_point, args=args)
                        elif obj == 1:
                            scipy_solve_point = self.secondary.polar_radius / 10.
                            solution = scipy.optimize.newton(self.secondary_potential_fn, scipy_solve_point, args=args)
                        # otestuje sa, ci vratena hodnota nie je NaN
                        if not np.isnan(solution):
                            # tu je taky logicky test, hviezda isto nemoze byt vacsi ako je vzdialenost medzi hviezdami
                            # a risenie urcite nemoze byt zaporne, kedze sa jedna o polomer vo sferickych suradniciach
                            # ak je vsetko yhovujuce, tak sa nastavi premenna pre pouzitie daneho riesenia na True
                            if actual_distance >= solution >= 0:
                                use = True
                    except:
                        use = False

                    # ak sa riesenie pouzije
                    if use:
                        # risenie sa ulozi do vektora a pretransformuje prislusnou funkciou zo sferickcyh do
                        # kartezianskych suradnic
                        point_coordinate[0], point_coordinate[1], point_coordinate[2] = solution, vector_spheric[1], \
                                                                                        vector_spheric[2]
                        xyz = Fn.spheric_to_cartesian(point_coordinate)
                        # ak sa jedna o primarnu zlozku, tak sa priamo pouzije x-ova suradnica
                        if obj == 0:
                            x = xyz[0]
                        # ak sa jedna o sekundarnu zlozku, je potrebne urobit urcite upravy, pretoze body sa
                        # pocitali po preneseni
                        # sekundarnej zlozky do pociatku suradnicoveho systemu, treba preto tuto zlozku posunut
                        # o vzdialenost hviezd
                        # a zrdkadlovo ju otocit
                        elif obj == 1:
                            x = - (xyz[0] - actual_distance)
                        y, z = xyz[1], xyz[2]
                        if transform == transform_steps - 1:
                            equatorial.append([x, y, z])
                        elif transform == 0:
                            meridional.append([x, y, z])
                        else:
                            partial.append([x, y, z])
                    # pricahdza na rad posunutie sa po kruznici v polomere "r"
                    # vektor sferickych suradnic, ktory ukazuje na dane miesto na kruznici sa pretrasformuje
                    # do kartezianskych suradnic
                    vector_xyz = Fn.spheric_to_cartesian(vector_spheric)
                    # nasledne sa karteziansky vektor ukazujuci na dane miesto na kruznici zrotuje o prislusncy
                    # uhol (zodpovedajuci
                    # dlzke kroku v danom polomere "r")
                    rotate = Fn.rotate(vector=vector_xyz, angle=rotation_angle)
                    # z kartezianskeho vektora uz prerotovaneho sa spatne ziska sfericky vektor, aby sme vedeli
                    # prislusny azimutalny
                    # a polarny uhol pre dalsie zratanie bodu na povrchu hviezdy
                    vector_spheric = Fn.cartesian_to_spheric(rotate)
                    # funkcia pre transformaciu moze vratit zaporny azimutalny uhol, tak sa to tu upravi,
                    # aby to bolo stale kladne
                    if vector_spheric[1] < 0: vector_spheric[1] += (2. * np.pi)

                # posun dalej v polarnom uhle
                theta += theta_step_length

            # funkcia ma na vstupe ovladaciu premennu, tzv. zero_point; tato premenna ovlada moznost vypoctu
            # bodu na povrchu
            # hviezdy v smere (phi = 0, theta = np.pi / 2.); je to takto riesene z dvoch dovodov
            # 1) pri krokovani v polarnom smere sa nemusi trafit funkcia presne na np.pi / 2.0
            # 2) ked sa spravne netrafi u napr. zlozky, ktora vyplna svoj rocheho lalok, tak model hviezdy
            # bude dost nekorektny
            # je to v podstate len prekopirovana cast vypoctu z kody vyssie
            if zero_point:
                # try count point with (r, 0, 180)
                args = (actual_distance, 0.0, np.pi / 2.)
                try:
                    # rozdiel oproti kodu vyssie je v tom, ze sa tu pouziva ina funkcia pre numericke riesenie
                    # korenov, pretoze fcia
                    # fsolve() disponuje obsiahlejsim vystupom a informaciami o vypocte, resp. an konci vypoctu,
                    # ktore je potreben overit
                    # snaha zratat tento bod pri wuam systemoch s tym, ze by sa tam naozaj nejaky bod nedopatrenim
                    # doratal (v hrdle),
                    # by viedla k nekorektnemu modelu
                    if obj == 0:
                        scipy_solve_point = self.primary.polar_radius / 10.
                        solution, info, ier, msg = scipy.optimize.fsolve(self.primary_potential_fn, scipy_solve_point,
                                                                         full_output=True, args=args)
                    elif obj == 1:
                        scipy_solve_point = self.secondary.polar_radius / 10.
                        solution, info, ier, msg = scipy.optimize.fsolve(self.secondary_potential_fn, scipy_solve_point,
                                                                         full_output=True, args=args)
                    if ier == 1 and not np.isnan(solution[0]):
                        solution = solution[0]
                        if actual_distance >= solution >= 0:
                            use = True
                except:
                    use = False
                if use:
                    point_coordinate[0], point_coordinate[1], point_coordinate[2] = solution, 0.0, np.pi / 2.
                    xyz = Fn.spheric_to_cartesian(point_coordinate)
                    if obj == 0:
                        x = xyz[0]
                    elif obj == 1:
                        x = - (xyz[0] - actual_distance)
                    y, z = xyz[1], xyz[2]
                    # riesenie sa ulozi do premennej z_point a nie priamo do pola partial, pretoze s tym sa potom bude
                    # este dalej pracovat tak, ze sa budu body preklapat a tento bod, by sa len duplikoval,
                    # co by pripadne
                    # pri pouziti cgal triangulacie sposobilo pad kodu
                    # partial.append([x, y, z])
                    z_point = [x, y, z]

            # prvy bod (phi = np.pi, theta = np.pi / 2.0) sa oddeli od rovnika, aby nedoslo k jeho zduplikovaniu (resp.
            # ku zrkadlovemu preklopeniu cez x os, co povedie k dvom bodom, ktore budu od seba vzdialene radovo 1e-15)
            f_point, equatorial = equatorial[0], equatorial[1:]

            # dopni sa prvy bod (first point), len raz, aby nedoslo k jeho zduplikovaniu
            # doplni sa zero point, aby nedoslo k jeho duplikacii

            full_partial = []
            if Fn.empty(z_point) and zero_point:
                print("ValueError: In class: Binary, function: get_3d_model_optimized(), line: " + str(
                    Fn.lineno()) + ". Variable `zero_point` is set to True, but there is not correct solution. "
                                   "Target object: " + str(obj) + ".")
                return False
            elif not Fn.empty(z_point) and zero_point:
                full_partial.append(z_point)

            if not Fn.empty(partial) and not Fn.empty(equatorial) and not Fn.empty(meridional):
                full_partial.append(f_point)
                for point in partial:
                    full_partial.append(point)
                    full_partial.append([point[0], -point[1], point[2]])
                    full_partial.append([point[0], point[1], -point[2]])
                    full_partial.append([point[0], -point[1], -point[2]])

                for point in equatorial:
                    full_partial.append(point)
                    full_partial.append([point[0], -point[1], point[2]])
                for point in meridional:
                    full_partial.append(point)
                    full_partial.append([point[0], point[1], -point[2]])

                total = total + full_partial

                if obj == 0:
                    model['primary'] = full_partial
                elif obj == 1:
                    model['secondary'] = full_partial
            else:
                print("ValueError: In class: Binary, function: get_3d_model_optimized(), line: " + str(
                    Fn.lineno()) + ". One of lists (`partial`, `equatorial`, `meridional`) is empty.")
                return False

        # computing of separation neck
        if self.binary_morph == "over-contact":
            separation, neck_radius = self.get_separation(actual_distance=actual_distance)
            r = min([self.primary.polar_radius, self.secondary.polar_radius])

            # cylindric phi angle steps and cylindric phi angle
            cps, cp = int(math.ceil(phi_steps * r / neck_radius)) * 2, 0.0
            step_length = 2.0 * np.pi / cps
            plane = []
            while True:
                if cp >= 2.0 * np.pi:
                    break
                args = actual_distance, cp, separation
                use = False
                solution = scipy.optimize.newton(self.wc_potential_fn,
                                                 min([self.primary.polar_radius / 10.0,
                                                      self.secondary.polar_radius / 10.0]),
                                                 args=args)
                # otestuje sa, ci vratena hodnota nie je NaN
                if not np.isnan(solution):
                    if max([self.primary.polar_radius, self.secondary.polar_radius]) >= solution >= 0:
                        use = True
                if use:
                    plane.append([separation, solution * np.cos(cp), solution * np.sin(cp)])
                cp += step_length

            # reorganise vertices (if any vertex of primary are located behind separation value put it to secondary
            # and vice-versa)
            primary, secondary = [], []
            for t in total:
                if t[0] > separation:
                    secondary.append(t)
                elif t[0] < separation:
                    primary.append(t)
            p, s = np.concatenate((primary, plane), 0), np.concatenate((secondary, plane), 0)
            u, v, w = len(p), len(s), len(plane)
            # average spacing
            avsp = [Fn.average_spacing(data=p, neighbours=4), Fn.average_spacing(data=s, neighbours=4)]

            remove, v, o, l = [[], []], [p, s], [None, None], [u, v]

            for t_o in [0, 1]:
                for plane_vertex in plane:
                    a = ([ndx for ndx, x in list(zip(range(0, l[t_o] - w), v[t_o][0:l[t_o] - w]))
                          if np.linalg.norm(np.array(plane_vertex) - np.array(x)) < avsp[t_o] / 2.0])
                    remove[t_o] = np.concatenate((remove[t_o], a), axis=0)

                remove[t_o] = np.array(list(set(remove[t_o])), dtype="int32").tolist()
                o[t_o] = np.array([x for ndx, x in list(zip(range(0, l[t_o] - w), v[t_o][0:l[t_o] - w]))
                                   if ndx not in remove[t_o]]).tolist()

            model["primary"] = np.concatenate((o[0], plane), 0)
            model["secondary"] = np.concatenate((o[1], plane), 0)
            model["system"] = np.concatenate((o[0], o[1], plane), 0)
            model["separation"] = separation, neck_radius

            # import Plot.Plot as Plt
            # Plt.plot_3d(vertices=[model["primary"]], faces=None, normals_view=False, points_view=True,
            #             faces_view=False, point_color="r", point_size=5.0, face_color="c", edge_color="k")

            del(primary, secondary, p, s, o, v, total)
            return model

        if total:
            model["separation"] = None
            model['system'] = np.array(total)
        return model

    def get_separation(self, actual_distance=None):
        linear_steps, x, phi = 0.005, 0.0, np.pi / 2.0
        separation = []
        while True:
            args = actual_distance, phi, x
            try:
                use = False
                solution = scipy.optimize.newton(self.wc_potential_fn,
                                                 min([self.primary.polar_radius / 10.0,
                                                      self.secondary.polar_radius / 10.0]),
                                                 args=args)
                # otestuje sa, ci vratena hodnota nie je NaN
                if not np.isnan(solution):
                    if max([self.primary.polar_radius, self.primary.polar_radius]) >= solution >= 0:
                        use = True
                if use:
                    separation.append([x, solution])
            except:
                pass

            x += linear_steps
            if x >= 1:
                break

        z = np.polyfit(list(zip(*separation))[0], list(zip(*separation))[1], 17)
        p = np.poly1d(z)
        # find minimum of function
        crit = p.deriv().r
        r_crit = crit[crit.imag == 0].real
        test = p.deriv(2)(r_crit)
        x_min = [d for d in r_crit[test > 0] if 0.0 < d < 1.0][0]

        # xs = np.arange(0., 1., 0.01)
        # ys = p(xs)

        # import matplotlib.pyplot as plt
        # plt.scatter(list(zip(*separation))[0], list(zip(*separation))[1], c="b")
        # plt.plot(xs, ys, c="r")
        # plt.axis('equal')
        # plt.show()

        return x_min, p(x_min)

class Orbit(object):
    def __init__(
            self,
            orbital_period=None,
            eccentricity=None,
            inclination=None,
            argument_of_periastron=None
    ):
        self.exception = []
        self.init = True
        self.inclination = None
        self.argument_of_periastron = argument_of_periastron

        if self.argument_of_periastron >= np.pi / 2.0:
            self.argument_of_periastron_azimut = argument_of_periastron - np.pi / 2.0
        else:
            self.argument_of_periastron_azimut = 2.0 * np.pi - self.argument_of_periastron

        if inclination < 0 or inclination > np.pi:
            self.exception.append("ValueError: In reference: Geometry, function: rotate_inclination, line: " + str(
                Fn.lineno()) + " Variable `inclination` is invalid. Use `inclination` in range [0, pi].")
            self.init = False
        else:
            self.inclination = inclination

        self.relative_semimajor_axis = None
        self.orbital_period = orbital_period
        self.eccentricity = eccentricity
        self.mean_angular_velocity = self.angular_velocity()  # $\rad.sec^{-1}$

        try:

            self.conjuction = self.conjunction(eccentricity=self.eccentricity,
                                               argument_of_periastron=self.argument_of_periastron)
        except (AttributeError, ValueError):
            self.init = False
            self.exception.append("ValueError: In class: Orbit, function: __init__(), line: " + str(
                Fn.lineno()) + ". Error has been occurred during initialisation, in function conjunction")

        # hodnoty pre primarny [0] a sekundarny [1] zakryt
        # prava anomalia konjunkcie (teda merana od periastra, od apsidalnej priamky)
        # \nu_{con}
        try:
            self.true_anomaly_of_conjuction = [self.conjuction[0]["true_anomaly"], self.conjuction[1]["true_anomaly"]]
        except (AttributeError, ValueError):
            self.init = False
            self.exception.append("ValueError: In class: Orbit, function: __init__(), line: " + str(
                Fn.lineno()) + ". Error has been occurred during initialisation, in function "
                               "true_anomaly_of_conjuction")

        # excentricka anomalia konjunkcie, merana od apsidalnej priamky
        try:
            self.eccentric_anomaly_of_conjunction = [self.conjuction[0]["eccentric_anomaly"],
                                                     self.conjuction[1]["eccentric_anomaly"]]
        except (AttributeError, ValueError):
            self.init = False
            self.exception.append("ValueError: In class: Orbit, function: __init__(), line: " + str(
                Fn.lineno()) + ". Error has been occurred during initialisation, in function "
                               "eccentric_anomaly_of_conjunction")

        # stredna anomalia konjunkcie, merana od apsidalnej priamky
        try:
            self.mean_anomaly_of_conjunction = [self.conjuction[0]["mean_anomaly"], self.conjuction[1]["mean_anomaly"]]
        except (AttributeError, ValueError):
            self.init = False
            self.exception.append("ValueError: In class: Orbit, function: __init__(), line: " + str(
                Fn.lineno()) + ". Error has been occurred during initialisation, in function "
                               "mean_anomaly_of_conjunction")

        # skutocna faza periastra a konjunkcie je totozna, neviem co som tym chcel povedat, ked som to programoval
        try:
            self.true_phase_of_conjunction = [self.conjuction[0]["true_phase"], self.conjuction[1]["true_phase"]]
        except (AttributeError, ValueError):
            self.init = False
            self.exception.append("ValueError: In class: Orbit, function: __init__(), line: " + str(
                Fn.lineno()) + ". Error has been occurred during initialisation, in function "
                               "true_phase_of_conjunction")

        # \Phi_{per}
        # faza periastra merana od nuly suradnioveho systemu po nulu sytemu orbitalnej drahy, v kladnom smere teda po
        # periastrum toto je vlastne shift na prevod medzi fazou orbitalnou a fotometrickou, kedze orbitalna je merana
        # od apsidalnej priamky
        try:
            self.true_phase_of_periastron = [self.conjuction[0]["true_phase"], self.conjuction[1]["true_phase"]]
        except (AttributeError, ValueError):
            self.init = False
            self.exception.append("ValueError: In class: Orbit, function: __init__(), line: " + str(
                Fn.lineno()) + ". Error has been occurred during initialisation, in function "
                               "true_phase_of_periastron")
        try:
            self.periastron_distance = self.orbital_motion(photometric_phase_from=0.0, photometric_phase_to=0.0,
                                                           shift=self.true_phase_of_periastron[0],
                                                           photometric_phase_step=1.0, eccentricity=self.eccentricity,
                                                           argument_of_periastron=self.argument_of_periastron)[0][0]
        except (AttributeError, ValueError):
            self.init = False
            self.exception.append("ValueError: In class: Orbit, function: __init__(), line: " + str(
                Fn.lineno()) + ". Error has been occurred during initialisation, in function "
                               "orbital_motion")

    def angular_velocity(self):
        return ((2.0 * np.pi) / (self.orbital_period * 86400.0)) * np.sqrt(
            (1.0 - self.eccentricity) * (1.0 + self.eccentricity))

    @classmethod
    def mean_anomaly(cls, phase=None):
        return 2.0 * np.pi * phase

    @classmethod
    def mean_anomaly_fn(cls, eccentric_anomaly, *args):
        mean_anomaly, eccentricity = args
        return (eccentric_anomaly - eccentricity * np.sin(eccentric_anomaly)) - mean_anomaly

    @classmethod
    def eccentric_anomaly(cls, eccentricity=None, mean_anomaly=None):
        import scipy.optimize
        args = (mean_anomaly, eccentricity)
        try:
            solution = scipy.optimize.newton(cls.mean_anomaly_fn, 1.0, args=args, tol=1e-10)

            if not np.isnan(solution):
                if solution < 0:
                    solution += 2.0 * np.pi
                return solution
            else:
                return False
        except:
            return False

    @classmethod
    def true_phase(cls, photometric_phase=None, shift=None):
        return photometric_phase + shift

    @classmethod
    def true_anomaly(cls, eccentricity=None, eccentric_anomaly=None):
        anomaly = 2.0 * np.arctan(
            np.sqrt((1.0 + eccentricity) / (1.0 - eccentricity)) * np.tan(eccentric_anomaly / 2.0))
        if anomaly < 0:
            anomaly += 2.0 * np.pi
        return anomaly

    @classmethod
    def radius(cls, true_anomaly=None, eccentricity=None):
        return (1.0 - eccentricity ** 2) / (1.0 + eccentricity * np.cos(true_anomaly))

    @classmethod
    def true_anomaly_to_azimuth(cls, true_anomaly=None, argument_of_periastron=None):
        azimut = np.float64(true_anomaly + argument_of_periastron - np.pi * 0.5)
        if azimut > (2.0 * np.pi):
            azimut -= 2.0 * np.pi
        if azimut < 0:
            azimut += 2.0 * np.pi
        return np.float64(azimut)

    @classmethod
    def conjunction(cls, argument_of_periastron=None, eccentricity=None):
        ret = {}
        for alpha, idx in list(zip([np.pi / 2.0, 3.0 * np.pi / 2.0], [0, 1])):
            # prava anomalia konjunkcie (teda merana od periastra, od apsidalnej priamky)

            true_anomaly_of_conjuction = alpha - argument_of_periastron  # \nu_{con}
            if true_anomaly_of_conjuction < 0:
                true_anomaly_of_conjuction += 2.0 * np.pi

            # excentricka anomalia konjunkcie, merana od apsidalnej priamky
            eccentric_anomaly_of_conjunction = 2.0 * np.arctan(
                np.sqrt((1.0 - eccentricity) / (1.0 + eccentricity)) * np.tan(true_anomaly_of_conjuction / 2.0))

            if eccentric_anomaly_of_conjunction < 0:
                eccentric_anomaly_of_conjunction += 2.0 * np.pi

            # stredna anomalia konjunkcie, merana od apsidalnej priamky
            mean_anomaly_of_conjunction = eccentric_anomaly_of_conjunction - eccentricity * np.sin(
                eccentric_anomaly_of_conjunction)
            if mean_anomaly_of_conjunction < 0:
                mean_anomaly_of_conjunction += 2.0 * np.pi

            true_phase_of_conjunction = mean_anomaly_of_conjunction / (2.0 * np.pi)

            if true_phase_of_conjunction < 0:
                true_phase_of_conjunction += 1.0

            ret[idx] = {}
            ret[idx]["true_anomaly"] = true_anomaly_of_conjuction
            ret[idx]["eccentric_anomaly"] = eccentric_anomaly_of_conjunction
            ret[idx]["mean_anomaly"] = mean_anomaly_of_conjunction
            ret[idx]["true_phase"] = true_phase_of_conjunction
        return ret

    @classmethod
    def orbital_motion(cls, photometric_phase_from=None, photometric_phase_to=None, shift=None,
                       photometric_phase_step=None, eccentricity=None, argument_of_periastron=None):

        phase, position = photometric_phase_from, []
        while True:
            if phase > photometric_phase_to: break

            true_phase = cls.true_phase(photometric_phase=phase, shift=shift)
            # toto prerobi cislo vacsie ako 1.0 na desatine, napr 1.1 na 0.1
            if abs(true_phase) > 1: true_phase = math.modf(true_phase)[0]
            if true_phase < 0: true_phase += 1.0

            mean_anomaly = cls.mean_anomaly(phase=true_phase)
            eccentric_anomaly = cls.eccentric_anomaly(eccentricity=eccentricity, mean_anomaly=mean_anomaly)
            true_anomaly = cls.true_anomaly(eccentric_anomaly=eccentric_anomaly, eccentricity=eccentricity)

            if true_anomaly < 0:
                true_anomaly = (2.0 * np.pi) - abs(true_anomaly)

            actual_distance = cls.radius(true_anomaly=true_anomaly, eccentricity=eccentricity)
            azimutal_angle = cls.true_anomaly_to_azimuth(true_anomaly=true_anomaly,
                                                         argument_of_periastron=argument_of_periastron)

            if azimutal_angle < 0:
                azimutal_angle = (2.0 * np.pi) - abs(azimutal_angle)

            position.append([actual_distance, azimutal_angle, true_anomaly, phase])

            phase += photometric_phase_step
        return position

    def get_periastron_distance(self):
        return self.periastron_distance

    def get_exception(self):
        return self.exception

    def get_eccentricity(self):
        return self.eccentricity


class Geometry(object):
    @staticmethod
    def normal_estimation(binary_object=None, actual_distance=None, vertices=None, t_object="primary",
                          mode="in_center"):
        warnings.simplefilter(action="ignore", category=FutureWarning)
        # kontrola premennych
        if binary_object is None:
            print("ValueError: In reference: Geometry, function: normal_estimation(), line: " + str(
                Fn.lineno()) + ". Binary object is set to None.")
            return False

        if Fn.numpy_array_is_empty_verbose(arr=vertices, function_name="normal_estimation", class_name="Geometry",
                                           var_name="vertices", verbose=True, line=str(Fn.lineno())):
            return False

        # ak je zle nastaveny objekt pre vypocet
        if t_object != "primary" and t_object != "secondary":
            print("ValueError: In reference: Geometry, function: normal_estimation(), line: " + str(
                Fn.lineno()) + ". Variable `object` is invalid, use `primary` or `secondary`.")
            return False

        # koniec kontroly premennych
        if t_object == 'primary':
            synchronicity_parameter = binary_object.primary.synchronicity_parameter
        elif t_object == 'secondary':
            synchronicity_parameter = binary_object.secondary.synchronicity_parameter
        mass_ratio, vector, index = binary_object.mass_ratio, [], 0

        # tie ifka su to preto, ze mam zadefinovany parameter mozno pre buduce pouzitie vo funkcii, a ked nie je
        # pouzity, tak ten priblby PyCharm pyskuje
        if mode == "in_point":
            pass
        if mode == "in_center":
            pass

        for point in vertices:
            try:
                x = point[0]
                r_square = x ** 2 + point[1] ** 2 + point[2] ** 2
                rw_square = (actual_distance - x) ** 2 + point[1] ** 2 + point[2] ** 2

                # vypocet derivacie podla premennej x
                # je potrebne rozlisovat primarnu a sekundarnu zlozku
                if t_object == "primary":
                    der_x = - (x / r_square ** (3.0 / 2.0)) + (
                        (mass_ratio * (actual_distance - x)) / rw_square ** (
                            3.0 / 2.0)) + synchronicity_parameter ** 2 * (
                        mass_ratio + 1) * x - (mass_ratio / actual_distance ** 2)
                else:
                    der_x = - (x / r_square ** (3.0 / 2.0)) + (
                        (mass_ratio * (actual_distance - x)) / rw_square ** (
                            3.0 / 2.0)) - synchronicity_parameter ** 2 * (
                        mass_ratio + 1) * (1. - x) + (1. / actual_distance ** 2)

                # pre derivacie podla premennej y a z ine je potrebne rozlisovat, ci sa pocitaju pre sekudnarnu alebo
                # pre primarnu zlozku
                der_y = - point[1] * ((1.0 / r_square ** (3.0 / 2.0)) + (mass_ratio / rw_square ** (3.0 / 2.0)) - (
                    (synchronicity_parameter ** 2) * (mass_ratio + 1)))
                der_z = - point[2] * ((1.0 / r_square ** (3.0 / 2.0)) + (mass_ratio / rw_square ** (3.0 / 2.0)))

                # vector.append([-der_x, -der_y, -der_z, index])
                vector.append([-der_x, -der_y, -der_z])

            except:
                print("Error: In reference: Geometry, function: normal_estimation(), line: " + str(
                    Fn.lineno()) + ". Cannot estimate normal vector in step " + str(index) + ".")
                return False
            index += 1

        return np.array(vector)

    @staticmethod
    def vector_array_translation(vector_arr=None, translation_arr=None):
        translated = [np.array(vector_arr[idx]) + np.array(translation_arr[idx]) for idx in range(len(vector_arr))]
        return np.array(translated)

    @staticmethod
    def face_orientation_a(face=None, t_object="primary", actual_distance=None):
        com, vec = Fn.center_of_mass(faces=face), []

        center = np.array([0., 0., 0.])
        if t_object == "secondary":
            center[0] = actual_distance

        for idx in range(0, len(com)):
            vec_a = (face[idx][0] - face[idx][1])
            vec_b = (face[idx][1] - face[idx][2])
            cross_p = np.cross(vec_a, vec_b)

            test_v = com[idx] - center
            if np.dot(test_v, cross_p) < 0:
                cross_p *= -1.0

            vec.append(cross_p)

        return np.array(vec)

    @staticmethod
    def angle(u=None, v=None):
        return np.arccos(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))

