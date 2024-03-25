import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.lang.Math;
import java.awt.Color;

public class Raytracer {
    // Define light direction
    private static final double LIGHT_X = 1.0;
    private static final double LIGHT_Y = 1.0;
    private static final double LIGHT_Z = 1.0;

    // Define light color
    private static final Color LIGHT_COLOR = Color.WHITE;

    // Define triangle vertices (Assuming counter-clockwise ordering)
    private static final Vector3D vertex1 = new Vector3D(0, 0, 0);
    private static final Vector3D vertex2 = new Vector3D(1, 0, 0);
    private static final Vector3D vertex3 = new Vector3D(0, 1, 0);

    // Define surface color
    private static final Color SURFACE_COLOR = Color.RED;

    // Image dimensions
    private static final int IMAGE_WIDTH = 500;
    private static final int IMAGE_HEIGHT = 500;

    public static void main(String[] args) throws IOException {
        String filename = "teapot.obj"; // Path to your OBJ file
        double[] point = {1, 2, 3}; // Example point
        double[] linePoint1 = {0, 0, 0}; // Example line point 1
        double[] linePoint2 = {1, 1, 1}; // Example line point 2
        double[] line1Point1 = {0, 0, 0};
        double[] line1Point2 = {1, 1, 1};
        double[] line2Point1 = {0, 0, 0};
        double[] line2Point2 = {1, -1, 1};

        // 1/2: Read an obj file and store it as a list
        try {
            List<float[]> vertices = readVertices(filename);
            System.out.println("Vertices:");
            for (float[] vertex : vertices) {
                System.out.println("(" + vertex[0] + ", " + vertex[1] + ", " + vertex[2] + ")");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // 3: Compute the distance between a point and a line
        double distance = distancePointLine(point, linePoint1, linePoint2);
        System.out.println("Distance between the point and the line: " + distance);

        // 4: Determine whether a line and triangle intersect


        // 5: Determine which faces from the obj file a line intersects


        // 6: Compute the normal of a triangle in 3D
        Vector3D normal = computeNormal(vertex1, vertex2, vertex3);
        System.out.println("Normal: " + normal); // Print the normal vector

        // 7: Compute the angle between two lines
        double angle = angleBetweenLines(line1Point1, line1Point2, line2Point1, line2Point2);
        System.out.println("Angle between the two lines: " + angle + " degrees");

        // 8: Compute the illumination for a triangle
        BufferedImage image = new BufferedImage(IMAGE_WIDTH, IMAGE_HEIGHT, BufferedImage.TYPE_INT_RGB);
        double intensity = computeLambertianIntensity(normal); // Compute intensity based on Lambertian shading
        Color illuminatedColor = applyLightIntensity(SURFACE_COLOR, intensity); // Apply light color and intensity to the surface color
        drawTriangle(image, vertex1, vertex2, vertex3, illuminatedColor); // Draw triangle with illuminated color

        // 9: Raytrace an obj in 3D for illumination intensities over a 2D screen


        // 10: Output screen pixels as an image
//        exportPNG();
    }

    //==================================================================================================================
    //TODO ========== 1/2: Read an obj file and store it as a list ==========
    public static List<float[]> readVertices(String filename) throws IOException {
        List<float[]> vertices = new ArrayList<>();
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line;

        while ((line = reader.readLine()) != null) {
            if (line.startsWith("v ")) {
                String[] parts = line.split("\\s+");
                float x = Float.parseFloat(parts[1]);
                float y = Float.parseFloat(parts[2]);
                float z = Float.parseFloat(parts[3]);
                vertices.add(new float[]{Math.abs(x), Math.abs(y), Math.abs(z)}); // added abs to get absolute value
            }
        }
        reader.close();
        return vertices;
    }

    //==================================================================================================================
    //TODO ========== 3: Compute the distance between a point and a line ==========
    public static double distancePointLine(double[] point, double[] linePoint1, double[] linePoint2) {
        double[] lineVector = new double[3];
        double[] pointToLineVector = new double[3];

        // Calculate the vector along the line
        for (int i = 0; i < 3; i++) {
            lineVector[i] = linePoint2[i] - linePoint1[i];
        }

        // Calculate the vector from a point on the line to the given point
        for (int i = 0; i < 3; i++) {
            pointToLineVector[i] = point[i] - linePoint1[i];
        }

        // Calculate the projection of the point-to-line vector onto the line vector
        double projection = dotProduct(pointToLineVector, lineVector) / dotProduct(lineVector, lineVector);

        // Calculate the closest point on the line to the given point
        double[] closestPoint = new double[3];
        for (int i = 0; i < 3; i++) {
            closestPoint[i] = linePoint1[i] + projection * lineVector[i];
        }

        // Calculate the distance between the given point and the closest point on the line
        return distance(point, closestPoint);
    }

    // Function to compute dot product of two vectors
    public static double dotProduct(double[] vector1, double[] vector2) {
        double result = 0;
        for (int i = 0; i < vector1.length; i++) {
            result += vector1[i] * vector2[i];
        }
        return result;
    }

    // Function to compute distance between two points
    public static double distance(double[] point1, double[] point2) {
        double sum = 0;
        for (int i = 0; i < point1.length; i++) {
            sum += Math.pow(point1[i] - point2[i], 2);
        }
        return Math.sqrt(sum);
    }

    //==================================================================================================================
    //TODO ========== 4: Determine whether a line and triangle intersect ==========


    //==================================================================================================================
    //TODO ========== 5: Determine which faces from the obj file a line intersects ==========


    //==================================================================================================================
    //TODO ========== 6: Compute the normal of a triangle in 3D ==========

    // Method to compute the normal vector of the triangle
    private static Vector3D computeNormal(Vector3D v1, Vector3D v2, Vector3D v3) {
        Vector3D edge1 = v2.subtract(v1);
        Vector3D edge2 = v3.subtract(v1);
        return edge1.crossProduct(edge2).normalize();
    }

    // Simple Vector3D class for vector operations
    private static class Vector3D {
        double x, y, z;
        public Vector3D(double x, double y, double z) {
            this.x = x;
            this.y = y;
            this.z = z;
        }
        public Vector3D subtract(Vector3D other) {
            return new Vector3D(this.x - other.x, this.y - other.y, this.z - other.z);
        }
        public Vector3D crossProduct(Vector3D other) {
            return new Vector3D(this.y * other.z - this.z * other.y,
                    this.z * other.x - this.x * other.z,
                    this.x * other.y - this.y * other.x);
        }
        public Vector3D normalize() {
            double magnitude = Math.sqrt(x * x + y * y + z * z);
            return new Vector3D(x / magnitude, y / magnitude, z / magnitude);
        }
        public double dotProduct(Vector3D other) {
            return this.x * other.x + this.y * other.y + this.z * other.z;
        }
        @Override
        public String toString() {
            return "(" + x + ", " + y + ", " + z + ")";
        }
    }

    //==================================================================================================================
    //TODO ========== 7: Compute the angle between two lines ==========
    public static double angleBetweenLines(double[] line1Point1, double[] line1Point2,
                                           double[] line2Point1, double[] line2Point2) {
        double[] line1Vector = new double[3];
        double[] line2Vector = new double[3];

        // Calculate vectors along the lines
        for (int i = 0; i < 3; i++) {
            line1Vector[i] = line1Point2[i] - line1Point1[i];
            line2Vector[i] = line2Point2[i] - line2Point1[i];
        }

        // Compute dot product of the two vectors
        double dotProduct = dotProduct(line1Vector, line2Vector);

        // Compute magnitudes of the vectors
        double magnitude1 = magnitude(line1Vector);
        double magnitude2 = magnitude(line2Vector);

        // Compute the angle in radians using the dot product and magnitudes
        double angle = Math.acos(dotProduct / (magnitude1 * magnitude2));

        // Convert angle from radians to degrees
        return Math.toDegrees(angle);
    }

    // Function to compute magnitude of a vector
    public static double magnitude(double[] vector) {
        double sumOfSquares = 0;
        for (double component : vector) {
            sumOfSquares += component * component;
        }
        return Math.sqrt(sumOfSquares);
    }

    //==================================================================================================================
    //TODO ========== 8: Compute the illumination for a triangle ==========

    // Method to compute Lambertian intensity
    private static double computeLambertianIntensity(Vector3D normal) {
        Vector3D lightDirection = new Vector3D(LIGHT_X, LIGHT_Y, LIGHT_Z).normalize();
        return Math.max(0, normal.dotProduct(lightDirection));
    }

    // Method to apply light intensity to surface color
    private static Color applyLightIntensity(Color surfaceColor, double intensity) {
        int red = (int) (surfaceColor.getRed() * intensity * LIGHT_COLOR.getRed() / 255);
        int green = (int) (surfaceColor.getGreen() * intensity * LIGHT_COLOR.getGreen() / 255);
        int blue = (int) (surfaceColor.getBlue() * intensity * LIGHT_COLOR.getBlue() / 255);
        return new Color(red, green, blue);
    }

    // Method to draw a triangle with specified color
    private static void drawTriangle(BufferedImage image, Vector3D v1, Vector3D v2, Vector3D v3, Color color) {
        // Implementation of triangle rasterization
        // (Not provided here)
    }

    //==================================================================================================================
    //TODO ========== 9: Raytrace an obj in 3D for illumination intensities over a 2D screen ==========


    //==================================================================================================================
    //TODO ========== 10: Output screen pixels as an image ==========
    public static void exportPNG() {
        // width of the image
        int width = 963;
        // height of the image
        int height = 640;
        // For storing image in RAM
        BufferedImage image = null;
        // READ IMAGE
        try {
            File input_file = new File("Goomba.png");
            // image file path create an object of
            // BufferedImage type and pass as parameter the
            // width,  height and image int
            // type. TYPE_INT_ARGB means that we are
            // representing the Alpha , Red, Green and Blue
            // component of the image pixel using 8 bit
            // integer value.
            image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
            // Reading input file
            image = ImageIO.read(input_file);
            System.out.println("Reading complete.");
        }
        catch (IOException e) {
            System.out.println("Error: " + e);
        }
        // WRITE IMAGE
        try {
            // Output file path
            File output_file = new File("output.png");
            // Writing to file taking type and path as
            ImageIO.write(image, "png", output_file);
            System.out.println("Writing complete.");
        }
        catch (IOException e) {
            System.out.println("Error: " + e);
        }
    }
}