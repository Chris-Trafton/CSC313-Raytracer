import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.lang.Math;
import java.util.Vector;

public class Raytracer {
    public static void main(String[] args) throws IOException {
        String filename = "teapot.obj"; // Path to your OBJ file
        double[] point = {1, 2, 3}; // Example point
        double[] linePoint1 = {0, 0, 0}; // Example line point 1
        double[] linePoint2 = {1, 1, 1}; // Example line point 2
        double[] line1Point1 = {0, 0, 0};
        double[] line1Point2 = {1, 1, 1};
        double[] line2Point1 = {0, 0, 0};
        double[] line2Point2 = {1, -1, 1};

        try {
            List<float[]> vertices = readVertices(filename);
            System.out.println("Vertices:");
            for (float[] vertex : vertices) {
                System.out.println("(" + vertex[0] + ", " + vertex[1] + ", " + vertex[2] + ")");
            }

            // Calculate distance between the given point and the line
            double distance = distancePointLine(point, linePoint1, linePoint2);
            System.out.println("Distance between the point and the line: " + distance);

            // Calculate the angle between two lines
            double angle = angleBetweenLines(line1Point1, line1Point2, line2Point1, line2Point2);
            System.out.println("Angle between the two lines: " + angle + " degrees");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    //==================================================================================================================
    //TODO ========== Reads an Obj File and Stores it as a List ==========
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
    //TODO ========== Computes the Distance Between a Point and a Line ==========
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
    //TODO ========== Determine whether a line and triangle intersect ==========
    public static boolean intersects(double[] P0, double[] P1, double[] V0, double[] V1, double[] V2) {
        double[] E1 = subtract(V1, V0);
        double[] E2 = subtract(V2, V0);
        double[] D = subtract(P1, P0);
        double[] P = crossProduct(D, E2);

        double det = dotProduct(E1, P);

        // check if the ray is parallel to the triangle
        if (det > -1e-6 && det < 1e-6) {
            return false;
        }

        double invDet = 1.0 / det;
        double[] T = subtract(P0, V0);
        double u = dotProduct(T, P) * invDet;

        // check if the intersection point is inside the triangle
        if (u < 0 || u > 1)
            return false;

        double[] Q = crossProduct(T, E1);
        double v = dotProduct(D, Q) * invDet;

        // check if the intersection point is inside the triangle
        if (v < 0 || u + v > 1)
            return false;

        double t = dotProduct(E2, Q) * invDet;

        // check if the intersection point is in front of the line's starting point
        if (t > 1e-6)
            return true;

        return false;
    }

    public static double[] subtract(double[] vector1, double[] vector2) {
        double[] ret = new double[vector1.length];
        for (int i = 0; i < vector1.length; i++) {
            double d = vector1[i] - vector2[i];
            ret[i] = d;
        }
//        double u = vector1[0] - vector2[0];
//        double v = vector1[1] - vector2[1];
//        double w = vector1[2] - vector2[2];
        return ret;
    }

    public static double[] crossProduct(double[] vector1, double[] vector2) {
        double[] ret = new double[vector1.length];
        ret[0] = (vector1[1] * vector2[2]) - (vector1[2] * vector2[1]);
        ret[1] = (vector1[2] * vector2[0]) - (vector1[0] * vector2[2]);
        ret[2] = (vector1[0] * vector2[1]) - (vector1[1] * vector2[0]);
        return ret;
    }

    //==================================================================================================================
    //TODO ========== Determine which faces from the obj file a line intersects ==========


    //==================================================================================================================
    //TODO ========== Computes the Angle Between Two Lines ==========
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
    //TODO ========== Compute the illumination for a triangle ==========


    //==================================================================================================================
    //TODO ========== Raytrace an obj in 3D for illumination intensities over a 2D screen ==========


    //==================================================================================================================
    //TODO ========== Output screen pixels as an image ==========
//    public static void exportPNG() {
//        // width of the image
//        int width = 963;
//        // height of the image
//        int height = 640;
//        // For storing image in RAM
//        BufferedImage image = null;
//        // READ IMAGE
//        try {
//            File input_file = new File("Goomba.png");
//            // image file path create an object of
//            // BufferedImage type and pass as parameter the
//            // width,  height and image int
//            // type. TYPE_INT_ARGB means that we are
//            // representing the Alpha , Red, Green and Blue
//            // component of the image pixel using 8 bit
//            // integer value.
//            image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
//            // Reading input file
//            image = ImageIO.read(input_file);
//            System.out.println("Reading complete.");
//        }
//        catch (IOException e) {
//            System.out.println("Error: " + e);
//        }
//        // WRITE IMAGE
//        try {
//            // Output file path
//            File output_file = new File("output.png");
//            // Writing to file taking type and path as
//            ImageIO.write(image, "png", output_file);
//            System.out.println("Writing complete.");
//        }
//        catch (IOException e) {
//            System.out.println("Error: " + e);
//        }
//    }
}