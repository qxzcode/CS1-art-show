"""
TODO: module/program docstring

Developed and tested with Python 3.7.
"""

import turtle, math, random
from typing import Tuple, Collection, Iterable, Hashable
Point3D = Tuple[float,float,float]
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
    """Linearly mix two colors. A mix_amount of 0.0 is color1, and 1.0 is color2."""
    return [(1-mix_amount)*v1 + mix_amount*v2 for v1, v2 in zip(color1, color2)]


# colors
UPPER_SKY_COLOR = (140/255, 188/255, 237/255) # also used for fog fading
LOWER_SKY_COLOR = (188/255, 220/255, 243/255)
SNOW_COLOR = (1.0, 1.0, 1.0)
TREE_COLOR = (97/255, 130/255, 101/255)
ROCK_COLOR = (0.51, 0.5, 0.57)

# other parameters
DEFAULT_FOG_FACTOR = 0.9
LIGHT_DIRECTION = normalize(3, 6, 1)
MIN_BRIGHTNESS = 0.3


class Camera:
    """A camera for drawing points in 3D space onto the 2D screen."""
    
    def __init__(self, pos: Point3D, theta_x: float, theta_y: float, theta_z: float,
                 zoom: float, fog_factor: float = DEFAULT_FOG_FACTOR):
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
        return self._scale * dx/dz, -self._scale * dy/dz, dz
    
    def draw_triangle(self, p1: Point3D, p2: Point3D, p3: Point3D, color: Color):
        """Draw a 2D triangle given its three vertices."""
        color_str = "#%02x%02x%02x" % tuple([round(255.0*x) for x in color])
        turtle.getcanvas().create_polygon(p1+p2+p3, fill=color_str)
    
    def compute_shaded_color(self, p1: Point3D, p2: Point3D, p3: Point3D, color: Color) -> Color:
        """Shade the color of a triangle according to a directional light."""
        # compute the normal vector
        ax, ay, az = p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2]
        bx, by, bz = p1[0]-p3[0], p1[1]-p3[1], p1[2]-p3[2]
        nx = ay*bz - az*by
        ny = az*bx - ax*bz
        nz = ax*by - ay*bx
        mag = math.sqrt(nx*nx + ny*ny + nz*nz)
        nx /= mag
        ny /= mag
        nz /= mag
        
        # compute the brightness multiplier
        lx, ly, lz = LIGHT_DIRECTION
        brightness = max(nx*lx + ny*ly + nz*lz, 0.0)
        brightness = MIN_BRIGHTNESS + brightness*(1-MIN_BRIGHTNESS)
        
        # return the darkened color
        return [v*brightness for v in color]
    
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
    """Handles generating and drawing fractal terrain."""
    
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
        fill_rectangle(min_x, min_y, max_x, min_y+y_step,
                       mix_colors(LOWER_SKY_COLOR, UPPER_SKY_COLOR, i/(num_steps-1)))
        min_y += y_step


def main():
    """The entry point of the program."""
    print("Setting up...")
    
    # set up turtle parameters
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
    fill_sky_gradient(50, 0.58)
    
    # set up the Camera
    camera = Camera((0, 6.0, -2.4), math.pi*0.34, 0, 0, zoom=3.4, fog_factor=0)
    camera = Camera((0, 0.07, -0.001), 0, 0, 0, zoom=1.2)
    
    # initialize the PRNG
    rand_seed = random.getrandbits(32)
    print("Seed:", rand_seed)
    
    # generate and draw the terrain
    print("Generating terrain...")
    random.seed(3911294863)
    color_offset = Terrain(recursion_depth=9, noise_depth=4, scale=0.35)
    random.seed(1462531765)
    terrain = Terrain(recursion_depth=9, noise_depth=7, scale=0.10,
                      snow_height=0.025, tree_height=-0.015, color_offset_heightmap=color_offset)
    terrain.draw(camera)
    print("Updating the screen...")
    turtle.update()
    
    # export the drawing to a file
    OUTPUT_FILE = "output.eps"
    print(f"Exporting {OUTPUT_FILE}...")
    turtle.getcanvas().postscript(file=OUTPUT_FILE, colormode="color", pagewidth="11i")
    
    # wait for the user to close the window
    print("Done!")
    turtle.mainloop()


if __name__ == "__main__":
    main()
