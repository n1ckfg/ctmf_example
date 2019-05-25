PImage depth, depth_filtered;

void setup() {
  size(1280, 480, P2D);
  depth = loadImage("shark.png");
  depth_filtered = depth.copy();
}

void draw() {
  image(depth, 0, 0);
  image(depth_filtered, 640, 0);
}
