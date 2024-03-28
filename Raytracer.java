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
        Vector3D point = new Vector3D(1, 2, 3);
        Vector3D linePoint1 = new Vector3D(0, 0, 0);
        Vector3D linePoint2 = new Vector3D(1, 1, 1);
        Vector3D line1Point1 = new Vector3D(0, 0, 0);
        Vector3D line1Point2 = new Vector3D(1, 1, 1);
        Vector3D line2Point1 = new Vector3D(0, 0, 0);
        Vector3D line2Point2 = new Vector3D(1, -1, 1);
        Vector3D lineStart = new Vector3D(0, 0, 0);
        Vector3D lineEnd = new Vector3D(1, 1, 1);

        // 1/2: Read an obj file and store it as a list
        List<Vector3D> vertices = new ArrayList<>();
        List<int[]> faces = new ArrayList<>();
        try {
            loadOBJ(filename, vertices, faces);
            System.out.println("Vertices:");
            for (Vector3D vertex : vertices) {
                System.out.println("(" + vertex.x + ", " + vertex.y + ", " + vertex.z + ")");
            }
            System.out.println("Faces:");
            for (int[] vertex : faces) {
                System.out.println("(" + vertex[0] + ", " + vertex[1] + ", " + vertex[2] + ")");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // 3: Compute the distance between a point and a line
        double distance = distancePointLine(point, linePoint1, linePoint2);
        System.out.println("Distance between the point and the line: " + distance);

        // 4: Determine whether a line and triangle intersect
        for (int i = 0; i < faces.size(); i += 3) {
            int[] face = faces.get(i);
            Vector3D V0 = vertices.get(face[0]); // vertices.get(i);
            Vector3D V1 = vertices.get(face[1]); // vertices.get(i + 1);
            Vector3D V2 = vertices.get(face[2]); // vertices.get(i + 2);

            // check for intersection between line and triangle
            boolean intersects = intersects(lineStart, lineEnd, V0, V1, V2);

            if (intersects) {
                System.out.println("Intersection detected with triangle " + (i / 3));
            }

            if( (i + 3) >= faces.size() - 3) {
                break;
            }
        }

        // 5: Determine which faces from the obj file a line intersects
        for (int i = 0; i < faces.size(); i++) {
            int[] face = faces.get(i);
            Vector3D v1 = vertices.get(face[0]);
            Vector3D v2 = vertices.get(face[1]);
            Vector3D v3 = vertices.get(face[2]);

            // Check intersection
            if (isIntersect(lineStart, lineEnd, v1, v2, v3)) {
                System.out.println("Line intersects with face " + i);
            }
        }

        // 6: Compute the normal of a triangle in 3D
        Vector3D normal = computeNormal(vertex1, vertex2, vertex3);
        System.out.println("Normal: " + normal); // Print the normal vector

        // 7: Compute the angle between two lines
        double angle = angleBetweenLines(line1Point1, line1Point2, line2Point1, line2Point2);
        System.out.println("Angle between the two lines: " + angle + " degrees");

        // 8: Compute the illumination for a triangle
        BufferedImage image = new BufferedImage(IMAGE_WIDTH, IMAGE_HEIGHT, BufferedImage.TYPE_INT_RGB); // TYPE_INT_ARGB
        double intensity = computeLambertianIntensity(normal); // Compute intensity based on Lambertian shading
        Color illuminatedColor = applyLightIntensity(SURFACE_COLOR, intensity); // Apply light color and intensity to the surface color
        drawTriangle(image, vertex1, vertex2, vertex3, illuminatedColor); // Draw triangle with illuminated color

        // 9: Raytrace an obj in 3D for illumination intensities over a 2D screen
        BufferedImage obj2D = raytraceObjTo2D(vertices);

        // 10: Output screen pixels as an image
        exportPNG(obj2D);
    }

    //==================================================================================================================
    //TODO ========== 1/2: Read an obj file and store it as a list ==========
    private static void loadOBJ(String filePath, List<Vector3D> vertices, List<int[]> faces) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(filePath));
        String line;
        while ((line = reader.readLine()) != null) {
            if (line.startsWith("v ")) {
                String[] parts = line.split("\\s+");
                double x = Double.parseDouble(parts[1]);
                double y = Double.parseDouble(parts[2]);
                double z = Double.parseDouble(parts[3]);
                vertices.add(new Vector3D(Math.abs(x), Math.abs(y), Math.abs(z)));
            } else if (line.startsWith("f ")) {
                String[] parts = line.split("\\s+");
                int[] face = new int[3];
                for (int i = 0; i < 3; i++) {
                    face[i] = Integer.parseInt(parts[i + 1]) - 1; // OBJ indices start from 1
                }
                faces.add(face);
            }
        }
        reader.close();
    }

    //==================================================================================================================
    //TODO ========== 3: Compute the distance between a point and a line ==========
    public static double distancePointLine(Vector3D point, Vector3D linePoint1, Vector3D linePoint2) {
        Vector3D lineVector = new Vector3D(0, 0, 0);
        Vector3D pointToLineVector = new Vector3D(0, 0, 0);

        // Calculate the vector along the line
        lineVector.x = linePoint2.x - linePoint1.x;
        lineVector.y = linePoint2.y - linePoint1.y;
        lineVector.z = linePoint2.z - linePoint1.z;

        // Calculate the vector from a point on the line to the given point
        pointToLineVector.x = point.x - linePoint1.x;
        pointToLineVector.y = point.y - linePoint1.y;
        pointToLineVector.z = point.z - linePoint1.z;

        // Calculate the projection of the point-to-line vector onto the line vector
        double projection = dotProduct(pointToLineVector, lineVector) / dotProduct(lineVector, lineVector);

        // Calculate the closest point on the line to the given point
        Vector3D closestPoint = new Vector3D(0, 0, 0);
        closestPoint.x = linePoint1.x + projection * lineVector.x;
        closestPoint.y = linePoint1.y + projection * lineVector.y;
        closestPoint.z = linePoint1.z + projection * lineVector.z;

        // Calculate the distance between the given point and the closest point on the line
        return distance(point, closestPoint);
    }

    // Function to compute dot product of two vectors
    public static double dotProduct(Vector3D vector1, Vector3D vector2) {
        double result = 0;
        result += vector1.x * vector2.x;
        result += vector1.y * vector2.y;
        result += vector1.z * vector2.z;

        return result;
    }

    // Function to compute distance between two points
    public static double distance(Vector3D point1, Vector3D point2) {
        double sum = 0;
        sum += Math.pow(point1.x - point2.x, 2);
        sum += Math.pow(point1.y - point2.y, 2);
        sum += Math.pow(point1.z - point2.z, 2);

        return Math.sqrt(sum);
    }

    //==================================================================================================================
    //TODO ========== 4: Determine whether a line and triangle intersect ==========

    // function given by ChatGPT and slightly modified for our use
    public static boolean intersects(Vector3D lineStart, Vector3D lineEnd, Vector3D V0, Vector3D V1, Vector3D V2) {
        Vector3D E1 = V1.subtract(V0);// subtract(V1, V0);
        Vector3D E2 = V2.subtract(V0); // subtract(V2, V0);
        Vector3D D = lineEnd.subtract(lineStart); // subtract(P1, P0);
        Vector3D P = D.crossProduct(E2); // crossProduct(D, E2);

        double det = E1.dotProduct(P); // dotProduct(E1, P);

        // check if the ray is parallel to the triangle
        if (det > -1e-6 && det < 1e-6) {
            return false;
        }

        double invDet = 1.0 / det;
        Vector3D T = lineStart.subtract(V0);// subtract(P0, V0);
        double u = T.dotProduct(P) * invDet; // dotProduct(T, P) * invDet;

        // check if the intersection point is inside the triangle
        if (u < 0 || u > 1)
            return false;

        Vector3D Q = T.crossProduct(E1); // crossProduct(T, E1);
        double v = D.dotProduct(Q) * invDet; // dotProduct(D, Q) * invDet;

        // check if the intersection point is inside the triangle
        if (v < 0 || u + v > 1)
            return false;

        double t = E2.dotProduct(Q) * invDet; // dotProduct(E2, Q) * invDet;

        // check if the intersection point is in front of the line's starting point
        if (t > 1e-6)
            return true;

        return false;
    }

    //==================================================================================================================
    //TODO ========== 5: Determine which faces from the obj file a line intersects ========== not complete
    public static boolean isIntersect(Vector3D lineStart, Vector3D lineEnd, Vector3D facePoint1, Vector3D facePoint2, Vector3D facePoint3) {
        // Calculate vectors
        double[] lineVector = {lineEnd.x - lineStart.x, lineEnd.y - lineStart.y, lineEnd.z - lineStart.z};

        // Check each edge of the triangle
        boolean intersects = intersect(lineStart, lineVector, facePoint1, facePoint2, facePoint3);
        intersects |= intersect(lineStart, lineVector, facePoint2, facePoint3, facePoint1);
        intersects |= intersect(lineStart, lineVector, facePoint3, facePoint1, facePoint2);

        return intersects;
    }

    // Check intersection between line and edge
    private static boolean intersect(Vector3D lineStart, double[] lineVector, Vector3D v1, Vector3D v2, Vector3D v3) {
        double[] edge1 = {v2.x - v1.x, v2.y - v1.y, v2.z - v1.z};
        double[] edge2 = {v3.x - v1.x, v3.y - v1.y, v3.z - v1.z};

        // Calculate normal to the triangle
        double[] normal = {
                edge1[1] * edge2[2] - edge1[2] * edge2[1],
                edge1[2] * edge2[0] - edge1[0] * edge2[2],
                edge1[0] * edge2[1] - edge1[1] * edge2[0]
        };

        // Calculate parameter t
        double dotProduct = normal[0] * lineVector[0] + normal[1] * lineVector[1] + normal[2] * lineVector[2];
        if (Math.abs(dotProduct) < 1e-8) {
            // Line is parallel to triangle, no intersection
            return false;
        }

        double[] startToV1 = {lineStart.x - v1.x, lineStart.y - v1.y, lineStart.z - v1.z};
        double t = -(normal[0] * startToV1[0] + normal[1] * startToV1[1] + normal[2] * startToV1[2]) / dotProduct;

        if (t < 0 || t > 1) {
            // Intersection point is outside line segment
            return false;
        }

        // Calculate intersection point
        double[] intersectionPoint = {
                lineStart.x + t * lineVector[0],
                lineStart.y + t * lineVector[1],
                lineStart.z + t * lineVector[2]
        };

        // Check if intersection point is inside the triangle
        double[] edge3 = {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
        double[] edge4 = {v2.x - v3.x, v2.y - v3.y, v2.z - v3.z};
        double[] edge5 = {v3.x - v1.x, v3.y - v1.y, v3.z - v1.z};

        double[] crossProduct1 = {
                edge1[1] * edge3[2] - edge1[2] * edge3[1],
                edge1[2] * edge3[0] - edge1[0] * edge3[2],
                edge1[0] * edge3[1] - edge1[1] * edge3[0]
        };
        double[] crossProduct2 = {
                edge2[1] * edge4[2] - edge2[2] * edge4[1],
                edge2[2] * edge4[0] - edge2[0] * edge4[2],
                edge2[0] * edge4[1] - edge2[1] * edge4[0]
        };
        double[] crossProduct3 = {
                edge3[1] * edge5[2] - edge3[2] * edge5[1],
                edge3[2] * edge5[0] - edge3[0] * edge5[2],
                edge3[0] * edge5[1] - edge3[1] * edge5[0]
        };

        double dotProduct1 = crossProduct1[0] * (intersectionPoint[0] - v1.x) + crossProduct1[1] * (intersectionPoint[1] - v1.y) + crossProduct1[2] * (intersectionPoint[2] - v1.z);
        double dotProduct2 = crossProduct2[0] * (intersectionPoint[0] - v2.x) + crossProduct2[1] * (intersectionPoint[1] - v2.y) + crossProduct2[2] * (intersectionPoint[2] - v2.z);
        double dotProduct3 = crossProduct3[0] * (intersectionPoint[0] - v3.x) + crossProduct3[1] * (intersectionPoint[1] - v3.y) + crossProduct3[2] * (intersectionPoint[2] - v3.z);

        return (dotProduct1 >= 0 && dotProduct2 >= 0 && dotProduct3 >= 0);
    }

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
    public static double angleBetweenLines(Vector3D line1Point1, Vector3D line1Point2, Vector3D line2Point1, Vector3D line2Point2) {
        Vector3D line1Vector = new Vector3D(0, 0, 0);
        Vector3D line2Vector = new Vector3D(0, 0, 0);

        // Calculate vectors along the lines
        line1Vector.x = line1Point2.x - line1Point1.x;
        line2Vector.x = line2Point2.x - line2Point1.x;
        line1Vector.y = line1Point2.y - line1Point1.y;
        line2Vector.y = line2Point2.y - line2Point1.y;
        line1Vector.z = line1Point2.z - line1Point1.z;
        line2Vector.z = line2Point2.z - line2Point1.z;

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
    public static double magnitude(Vector3D vector) {
        double sumOfSquares = 0;
        double[] tempVector = {vector.x, vector.y, vector.z};
        for (double component : tempVector) {
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
    public static BufferedImage raytraceObjTo2D(List<Vector3D> vertices) {
        BufferedImage ret = new BufferedImage(IMAGE_WIDTH, IMAGE_HEIGHT, BufferedImage.TYPE_INT_RGB);
        // first set all pixels to white
        for (int i = 0; i < IMAGE_WIDTH; i++) {
            for (int j = 0; j < IMAGE_HEIGHT; j++) {
                ret.setRGB(i, j, 0);
            }
        }

        // where a vertex is, set that pixel to red
        for (int i = 0; i < vertices.size(); i++) {
            Vector3D vertex = vertices.get(i);
            double x = vertex.x;
            double y = vertex.y;
            String color = SURFACE_COLOR.toString();
            int equalsIndex = color.indexOf("=") + 1;
            String rgb = color.substring(equalsIndex, equalsIndex + 3);
            int colorToInt = 65536 * Integer.parseInt(rgb) + 256 * 0 + 0; // Integer.parseInt(rgb);
            ret.setRGB((int)x, (int) y, colorToInt);
        }

        return ret;
    }

    //==================================================================================================================
    //TODO ========== 10: Output screen pixels as an image ==========
    public static void exportPNG(BufferedImage input) {
        // width of the image
        int width = IMAGE_WIDTH;
        // height of the image
        int height = IMAGE_HEIGHT;
        System.out.println("Reading complete.");
        // WRITE IMAGE
        try {
            // Output file path
            File output_file = new File("output.png");
            // Writing to file taking type and path as
            ImageIO.write(input, "png", output_file);
            System.out.println("Writing complete.");
        }
        catch (IOException e) {
            System.out.println("Error: " + e);
        }
    }

}