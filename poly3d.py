"""
Art show submission - draws a mountain range.

Developed and tested with Python 3.7.

Author: Quinn Tucker
"""

import turtle, math, random, argparse
from typing import Tuple, Collection, Iterable, Hashable
Point3D = Tuple[float,float,float]
Point2D = Tuple[float,float]
Color = Tuple[float,float,float]


def normalize(x: float, y: float, z: float) -> Point3D:
    """Normalize the given vector to be of unit length."""
    mag = math.sqrt(x*x + y*y + z*z)
    return x/mag, y/mag, z/mag


def progress_iterator(collection: Collection, message: str) -> Iterable:
    """Wrap a collection to print iteration progress as a percentage."""
    num_items = len(collection)
    last_percentage = -1
    for i, item in enumerate(collection):
        percentage = 100 * i // num_items
        if percentage > last_percentage:
            last_percentage = percentage
            print(f"{message} {percentage}%", end='\r')
        yield item
    print(f"{message} 100%")


def mix_colors(color1: Color, color2: Color, mix_amount: float) -> Color:
    """Linearly mix two colors. A mix_amount of 0.0 gives color1, and 1.0 gives color2."""
    return [(1-mix_amount)*v1 + mix_amount*v2 for v1, v2 in zip(color1, color2)]


def multiply_vectors(vec1: Iterable[float], vec2: Iterable[float]) -> Iterable[float]:
    """Multiply two vectors element-wise."""
    return [v1*v2 for v1, v2 in zip(vec1, vec2)]


# colors
UPPER_SKY_COLOR = (91/255, 130/255, 185/255) # also used for fog fading
LOWER_SKY_COLOR = (238/255, 207/255, 94/255)
SNOW_COLOR = (1.0, 1.0, 1.0)
TREE_COLOR = (97/255, 130/255, 101/255)
ROCK_COLOR = (49/94, 49/119, 61/143)
SUNLIGHT_COLOR = (255/255*1.5, 143/255*1.5, 86/255*1.5)
AMBIENT_LIGHT_COLOR = (87/255*0.2, 72/255*0.2, 93/255*0.2) # the color of shadows

# other parameters
FOG_FACTOR = 0.9
SUNLIGHT_DIRECTION = (-5, 1, 17)
AMBIENT_LIGHT_DIRECTION = (0, 1, 0)


class DirectionalLight:
    """A colored light that hits all surfaces from a constant direction."""
    
    def __init__(self, direction: Point3D, color: Color, dot_clip: float = 0.0):
        """
        Create a directional light given its direction and color.
        
        The dot_clip parameter adjusts the value of the dot product used
        in the lighting calculation; a lower value compresses the range of
        brightnesses produced by the light.
        """
        self._direction = normalize(*direction)
        self._color = color
        self._dot_clip = dot_clip
    
    def get_max_brightness(self) -> float:
        """Return the maximum color value that this light can produce."""
        return max(self._color)
    
    def compute_shaded_color(self, normal: Point3D, material_color: Color) -> Color:
        """
        Return the color contributed by this light on a surface given
        its (unit) normal vector and material color.
        """
        dot_product = sum(multiply_vectors(self._direction, normal))
        light_amount = max(dot_product, self._dot_clip)
        light_amount = (light_amount - self._dot_clip) / (1.0 - self._dot_clip)
        return [vm*vl*light_amount for vm, vl in zip(material_color, self._color)]


class Camera:
    """A camera for drawing points in 3D space onto the 2D screen."""
    
    def __init__(self, pos: Point3D, theta_x: float, theta_y: float, theta_z: float,
                 zoom: float, fog_factor: float,
                 lights: Iterable[DirectionalLight], fast_draw: bool):
        """Create a camera given its position in 3D space and the Tait-Bryan angles of its orientation."""
        self._pos = pos
        self._cx = math.cos(theta_x)
        self._sx = math.sin(theta_x)
        self._cy = math.cos(theta_y)
        self._sy = math.sin(theta_y)
        self._cz = math.cos(theta_z)
        self._sz = math.sin(theta_z)
        self._scale = turtle.window_height() * zoom
        self._fog_factor = fog_factor
        self._lights = lights
        self._fast_draw = fast_draw
    
    def project_point(self, point: Point3D) -> Point3D:
        """Project a point in 3D world space into 2D screen space."""
        x, y, z = point
        cam_x, cam_y, cam_z = self._pos
        x -= cam_x
        y -= cam_y
        z -= cam_z
        dx = self._cy*(self._sz*y + self._cz*x) - self._sy*z
        dy = self._sx*(self._sy*(self._sz*y + self._cz*x) + self._cy*z) + self._cx*(self._cz*y - self._sz*x)
        dz = self._cx*(self._sy*(self._sz*y + self._cz*x) + self._cy*z) - self._sx*(self._cz*y - self._sz*x)
        return self._scale * dx/dz, self._scale * dy/dz, dz
    
    def draw_triangle(self, p1: Point2D, p2: Point2D, p3: Point2D, color: Color):
        """Draw a 2D triangle given its three vertices."""
        if self._fast_draw:
            color_str = "#%02x%02x%02x" % tuple([round(255.0*x) for x in color])
            x1, y1 = p1
            x2, y2 = p2
            x3, y3 = p3
            turtle.getcanvas().create_polygon((x1,-y1,x2,-y2,x3,-y3), fill=color_str)
        else:
            turtle.goto(*p1)
            turtle.fillcolor(color)
            turtle.begin_fill()
            turtle.goto(*p2)
            turtle.goto(*p3)
            turtle.end_fill()
    
    def compute_shaded_color(self, p1: Point3D, p2: Point3D, p3: Point3D, material_color: Color) -> Color:
        """Shade the color of a triangle according to a directional light."""
        # compute the normal vector
        ax, ay, az = p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2]
        bx, by, bz = p1[0]-p3[0], p1[1]-p3[1], p1[2]-p3[2]
        nx = ay*bz - az*by
        ny = az*bx - ax*bz
        nz = ax*by - ay*bx
        normal = normalize(nx, ny, nz)
        
        # compute the color of the lights on the surface
        cr, cg, cb = 0, 0, 0
        for light in self._lights:
            lr, lg, lb = light.compute_shaded_color(normal, material_color)
            cr += lr
            cg += lg
            cb += lb
        
        # clip the color values at 1.0 and return
        max_v = max(cr, cg, cb)
        if max_v > 1.0:
            return min(cr,1.0), min(cg,1.0), min(cb,1.0)
        return cr, cg, cb
    
    def compute_fog_faded_color(self, color: Color, dz: float) -> Color:
        """Fade a color depending on how far from the camera it is."""
        fade_amount = math.exp(-(dz * self._fog_factor)**2)
        return mix_colors(UPPER_SKY_COLOR, color, fade_amount)
    
    def draw_triangles(self, triangles: Collection):
        """Shade, project, and draw a list of triangles in 3D."""
        # project the points into 2D and compute each shaded/faded color
        processed_triangles = []
        for p1, p2, p3, color in progress_iterator(triangles, "Processing triangles..."):
            shaded_color = self.compute_shaded_color(p1, p2, p3, color)
            *p1_p, z1 = self.project_point(p1)
            *p2_p, z2 = self.project_point(p2)
            *p3_p, z3 = self.project_point(p3)
            centroid_z = (z1 + z2 + z3) / 3
            faded_color = self.compute_fog_faded_color(shaded_color, centroid_z)
            processed_triangles.append((centroid_z, p1_p, p2_p, p3_p, faded_color))
        
        # sort the list of triangles back-to-front (by centroid Z depth)
        processed_triangles.sort(key=lambda tri: tri[0], reverse=True)
        
        # draw the triangles
        for _, p1, p2, p3, color in progress_iterator(processed_triangles, "Adding triangles to the canvas..."):
            self.draw_triangle(p1, p2, p3, color)
        print(f"    Added {len(processed_triangles)} triangles")


class Terrain:
    """Handles generating and drawing fractal heightmaps/terrain."""
    
    def __init__(self, recursion_depth: int, noise_depth: int, scale: float,
                 snow_height: float = None, tree_height: float = None,
                 color_offset_heightmap: "Terrain" = None):
        """Initialize and generate the Terrain."""
        self._depth = recursion_depth
        self._noise_depth = noise_depth
        self._scale = scale
        self._snow_height = snow_height
        self._tree_height = tree_height
        self._color_offset_heightmap = color_offset_heightmap
        self._height_cache = {}
        
        self._only_heightmap = snow_height is None or tree_height is None
        if self._only_heightmap:
            self._heightmap = {}
        else:
            self._triangles = []
        
        # generate the terrain recursively
        self._fractal_triangle((-1,0,math.sqrt(3)), (1,0,math.sqrt(3)), (0,0,0), depth=self._depth)
    
    def _get_midpoint(self, p1: Point3D, p2: Point3D, displace: bool) -> Point3D:
        """Return the two given points' midpoint, optionally displacing it randomly."""
        key = (p1, p2) if p1 < p2 else (p2, p1)
        if key not in self._height_cache:
            x1, y1, z1 = p1
            x2, y2, z2 = p2
            if displace:
                displacement = random.gauss(0, math.hypot(x1-x2, z1-z2) * self._scale)
            else:
                displacement = 0
            self._height_cache[key] = (x1+x2)/2, (y1+y2)/2 + displacement, (z1+z2)/2
        return self._height_cache[key]
    
    def _get_heightmap_key(self, p1: Point3D, p2: Point3D, p3: Point3D) -> Hashable:
        """Return the key in the _heightmap dict for the given triangle."""
        return p1[0]+p2[0]+p3[0], p1[2]+p2[2]+p3[2]
    
    def _fractal_triangle(self, p1: Point3D, p2: Point3D, p3: Point3D, depth: int):
        """Recursively subdivide the triangle, building the triangle list."""
        if depth == 0:
            height = (p1[1]+p2[1]+p3[1])/3
            if self._only_heightmap:
                self._heightmap[self._get_heightmap_key(p1,p2,p3)] = height
            else:
                if self._color_offset_heightmap is not None:
                    height += self._color_offset_heightmap.get_height(p1, p2, p3)
                if height > self._snow_height:
                    c = SNOW_COLOR
                elif height < self._tree_height:
                    c = TREE_COLOR
                else:
                    c = ROCK_COLOR
                self._triangles.append((p1, p2, p3, c))
        else:
            displace = depth <= self._noise_depth
            mid12 = self._get_midpoint(p1, p2, displace)
            mid23 = self._get_midpoint(p2, p3, displace)
            mid13 = self._get_midpoint(p3, p1, displace)
            self._fractal_triangle(p1, mid12, mid13, depth=depth-1)
            self._fractal_triangle(mid12, p2, mid23, depth=depth-1)
            self._fractal_triangle(mid13, mid23, p3, depth=depth-1)
            self._fractal_triangle(mid12, mid23, mid13, depth=depth-1)
    
    def get_height(self, p1: Point3D, p2: Point3D, p3: Point3D) -> float:
        """Return the height value for the given triangle."""
        return self._heightmap[self._get_heightmap_key(p1,p2,p3)]
    
    def draw(self, camera: Camera):
        """Draw the terrain using the given Camera."""
        camera.draw_triangles(self._triangles)


def fill_rectangle(min_x: float, min_y: float, max_x: float, max_y: float, color: Color):
    """Fill the axis-aligned rectangle bounded by the given coordinates."""
    turtle.goto(min_x, min_y)
    turtle.fillcolor(color)
    turtle.begin_fill()
    turtle.goto(max_x, min_y)
    turtle.goto(max_x, max_y)
    turtle.goto(min_x, max_y)
    turtle.end_fill()


def fill_sky_gradient(num_steps: int, start_y: float):
    """
    Fill the background sky gradient.
    
    Uses num_steps rectangles to approximate a linear gradient
    that goes from the top of the screen to start_y of the way
    down the screen (between 0.0 and 1.0).
    """
    # compute some helper values
    min_x = -turtle.window_width() / 2
    max_x = +turtle.window_width() / 2
    y_step = turtle.window_height()*start_y / num_steps
    min_y = turtle.window_height() / 2 - turtle.window_height()*start_y
    
    # fill the section below the gradient
    fill_rectangle(min_x, -turtle.window_height()/2, max_x, min_y, LOWER_SKY_COLOR)
    
    # fill the gradient
    for i in range(num_steps):
        fill_rectangle(min_x, min_y, max_x, min_y + y_step + 1,
                       mix_colors(LOWER_SKY_COLOR, UPPER_SKY_COLOR, i/(num_steps-1)))
        min_y += y_step


def main():
    """The entry point of the program."""
    # parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-export", action="store_true",
                        help="Don't export an .eps file of the drawing")
    parser.add_argument("--fast", action="store_true",
                        help="Add triangles directly to the Tkinter canvas for speed")
    parser.add_argument("--birds-eye", action="store_true",
                        help="Show a bird's eye view of the entire terrain")
    parser.add_argument("--random-terrain", action="store_true",
                        help="Use a random seed for the terrain heightmap")
    parser.add_argument("--random-color-offset", action="store_true",
                        help="Use a random seed for the color offset heightmap")
    args = parser.parse_args()
    
    # set up turtle parameters
    print("Setting up...")
    turtle.setup(9999, 9999)
    win_scale = min(turtle.window_width()//22, turtle.window_height()//17)
    turtle.setup(win_scale*22, win_scale*17) # the largest 11x8.5 window possible
    turtle.title("Submission by Quinn Tucker")
    turtle.tracer(0, 0)
    turtle.setundobuffer(None)
    turtle.hideturtle()
    turtle.penup()
    
    # fill the background with the sky gradient
    print("Filling the sky...")
    fill_sky_gradient(256, 0.58)
    
    # set up the lights and camera
    lights = [
        #DirectionalLight(SUNLIGHT_DIRECTION, SUNLIGHT_COLOR, dot_clip=0.0),
        DirectionalLight(AMBIENT_LIGHT_DIRECTION, AMBIENT_LIGHT_COLOR, dot_clip=-0.0),
    ]
    if args.birds_eye:
        camera = Camera((0, 6.0, -2.4), math.pi*0.34, 0, 0, zoom=3.4, fog_factor=0, lights=lights, fast_draw=args.fast)
    else:
        camera = Camera((0, 0.07, -0.001), 0, 0, 0, zoom=1.2, fog_factor=FOG_FACTOR, lights=lights, fast_draw=args.fast)
    
    # generate and draw the terrain
    print("Generating terrain...")
    if args.random_color_offset:
        color_offset_seed = random.getrandbits(32)
        print(f"    Color offset seed = {color_offset_seed}")
    else:
        color_offset_seed = 3038607546
    random.seed(color_offset_seed)
    color_offset = Terrain(recursion_depth=9, noise_depth=4, scale=0.35)
    
    if args.random_terrain:
        terrain_seed = random.getrandbits(32)
        print(f"    Terrain seed = {terrain_seed}")
    else:
        terrain_seed = 129477298
    random.seed(terrain_seed)
    terrain = Terrain(recursion_depth=9, noise_depth=7, scale=0.10,
                      snow_height=0.025, tree_height=-0.015, color_offset_heightmap=color_offset)
    
    terrain.draw(camera)
    print("Updating the screen...")
    turtle.update()
    
    # export the drawing to a file
    if not args.no_export:
        OUTPUT_FILE = "output.eps"
        print(f"Exporting {OUTPUT_FILE}...")
        turtle.getcanvas().postscript(file=OUTPUT_FILE, colormode="color", pagewidth="11i")
    
    # wait for the user to close the window
    print("Done!")
    turtle.mainloop()


if __name__ == "__main__":
    main()
