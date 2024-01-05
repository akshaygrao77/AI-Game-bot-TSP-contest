package root;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

public class MainDriver {
	public static void main(String[] args) throws IOException {
		BufferedReader bi = new BufferedReader(new InputStreamReader(System.in));

		String euclOrNonEucl = bi.readLine().trim();
		boolean isEuclideanDistance = euclOrNonEucl.equals("euclidean");

		int numberOfCities = Integer.parseInt(bi.readLine().trim());

		Node[] nodeList = new Node[numberOfCities];
		for (int i = 0; i < numberOfCities; i++) {
			String coordinates[] = bi.readLine().trim().split(" ");
			Node newNode = new Node();
			newNode.setX(Float.parseFloat(coordinates[0]));
			newNode.setY(Float.parseFloat(coordinates[1]));
			nodeList[i] = newNode;
		}

		float[][] distanceMatrix = new float[numberOfCities][numberOfCities];
		float[][] euclideanDistanceMatrix = new float[numberOfCities][numberOfCities];
		for (int i = 0; i < numberOfCities; i++) {
			String distancePerCity[] = bi.readLine().trim().split("\\s+");
			for (int j = 0; j < numberOfCities; j++) {
				distanceMatrix[i][j] = Float.parseFloat(distancePerCity[j]);
				if (isEuclideanDistance) {
					euclideanDistanceMatrix[i][j] = distanceMatrix[i][j];
				} else {
					euclideanDistanceMatrix[i][j] = calculateEuclideanDistance(nodeList[i], nodeList[j]);
				}
			}
		}
		TravellingSalesmanSolver tspSolver = new TravellingSalesmanSolver();
		tspSolver.solveProblem(numberOfCities, nodeList, isEuclideanDistance, distanceMatrix, euclideanDistanceMatrix);
	}

	private static float calculateEuclideanDistance(Node node1, Node node2) {
		float xDiff = Math.abs(node1.getX() - node2.getX());
		float xDiffSq = xDiff * xDiff;
		float yDiff = Math.abs(node1.getY() - node2.getY());
		float yDiffSq = yDiff * yDiff;

		return (float) Math.sqrt(xDiffSq + yDiffSq);
	}
}