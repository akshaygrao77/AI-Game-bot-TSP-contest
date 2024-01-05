package root;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Random;
import java.util.Set;

public class TravellingSalesmanSolver {
	public float leastTourCost = Float.MAX_VALUE;
	public ArrayList<Integer> leastCostTour = null;

	public void solveProblem(int numberOfCities, Node[] nodeList, boolean isEuclideanDistance, float[][] distanceMatrix,
			float[][] euclideanDistanceMatrix) {
		ArrayList<Integer> greedyTour = greedyHeuristics(isEuclideanDistance, numberOfCities, euclideanDistanceMatrix,
				distanceMatrix);

		ArrayList<ArrayList<Integer>> listOfToursByNearestNeighbours = nearestNeighbour(isEuclideanDistance,
				numberOfCities, euclideanDistanceMatrix, distanceMatrix);

		ArrayList<ArrayList<Integer>> listOfToursBySavingsHeuristics = savingsHeuristics(isEuclideanDistance,
				numberOfCities, euclideanDistanceMatrix, distanceMatrix);

		ArrayList<ArrayList<Integer>> accumulatedOverallTours = new ArrayList<ArrayList<Integer>>();
		accumulatedOverallTours.addAll(listOfToursByNearestNeighbours);
		accumulatedOverallTours.addAll(listOfToursBySavingsHeuristics);
		accumulatedOverallTours.add(greedyTour);
		int populationSize = numberOfCities * 5;
		while (true) {

			geneticAlgorithm(accumulatedOverallTours, numberOfCities, isEuclideanDistance, distanceMatrix,
					euclideanDistanceMatrix, populationSize);
			populationSize = populationSize * 2;
		}
	}

	private ArrayList<Integer> realignTourToCommonBase(ArrayList<Integer> tour) {
		int basePoint = 0;
		LinkedList<Integer> realignedTour = new LinkedList<Integer>();
		if (tour.get(0) == basePoint) {
			return new ArrayList<Integer>(tour);
		}
		int i = 0;
		while (tour.get(i) != basePoint) {
			realignedTour.addLast(tour.get(i));
			i++;
		}
		int j = tour.size() - 1;
		while (tour.get(j) != basePoint) {
			realignedTour.addFirst(tour.get(j));
			j--;
		}
		realignedTour.addFirst(basePoint);

		return new ArrayList<Integer>(realignedTour);
	}

	private void geneticAlgorithm(ArrayList<ArrayList<Integer>> accumulatedOverallTours, int numberOfCities,
			boolean isEuclideanDistance, float[][] distanceMatrix, float[][] euclideanDistanceMatrix,
			int initialPopulationSize) {
		int numberOfMutations = numberOfCities;
		int sizeOfToBeGenaratedPopulation = initialPopulationSize - accumulatedOverallTours.size();
		ArrayList<ArrayList<Integer>> initialPopulation = generateInitialPopulationOfTours(
				sizeOfToBeGenaratedPopulation, numberOfCities, accumulatedOverallTours);
		ArrayList<ArrayList<Integer>> currentPopulation = initialPopulation;
		ArrayList<Tour> listOfFitnessCalculatedTours = new ArrayList<Tour>(currentPopulation.size());
		for (ArrayList<Integer> eachTour : currentPopulation) {
			Tour newTour = new Tour();
			newTour.pathOfTour = eachTour;
			listOfFitnessCalculatedTours.add(newTour);
		}
		float leastPathCost = injectFitnessScoresToTour(distanceMatrix, currentPopulation.size(),
				listOfFitnessCalculatedTours);
		float previousPathCost = leastPathCost;
		float leastPathCostOffSpring = leastPathCost;
		for (int i = 0; i < numberOfMutations; i++) {
			do {
				previousPathCost = leastPathCostOffSpring;
				int selectedPopulationSize = listOfFitnessCalculatedTours.size() / 2;
				Collections.sort(listOfFitnessCalculatedTours);
				ArrayList<Tour> populationSelectedBasedOnFitness = selectPopulationBasedOnFitnessAndRandomness(
						selectedPopulationSize, listOfFitnessCalculatedTours);
				Collections.sort(populationSelectedBasedOnFitness);
				ArrayList<Tour> offSpringPopulation = crossOverSelectedPopulation(numberOfCities,
						populationSelectedBasedOnFitness);
				leastPathCostOffSpring = injectFitnessScoresToTour(distanceMatrix, offSpringPopulation.size(),
						offSpringPopulation);
				Collections.sort(offSpringPopulation);
				listOfFitnessCalculatedTours = mergeOffSpringAndOriginalPopulation(numberOfCities, offSpringPopulation,
						listOfFitnessCalculatedTours);
			} while (previousPathCost != leastPathCostOffSpring);
			// Initially mutation would be 25% of population.
			double mutationRatio = 0.25;
			if (i > (numberOfMutations / 2)) {
				mutationRatio = 0.5;
			}
			int numberOfToursToMutate = (int) (listOfFitnessCalculatedTours.size() * mutationRatio);
			randomlyMutateThePopulation(listOfFitnessCalculatedTours, numberOfToursToMutate, distanceMatrix);
		}
	}

	private void randomlyMutateThePopulation(ArrayList<Tour> listOfFitnessCalculatedTours, int numberOfToursToBeMutated,
			float[][] distanceMatrix) {
		HashSet<Tour> populationSet = new HashSet<Tour>(listOfFitnessCalculatedTours);
		HashSet<Integer> randomIndexSet = new HashSet<Integer>(numberOfToursToBeMutated);
		ArrayList<Integer> randomTour = listOfFitnessCalculatedTours.get(0).pathOfTour;
		int addCount = 0;
		while (addCount != numberOfToursToBeMutated) {
			Collections.shuffle(randomTour);
			Tour randomTourObj = new Tour();
			randomTourObj.pathOfTour = realignTourToCommonBase((ArrayList<Integer>) randomTour.clone());
			randomTourObj.pathCost = calculateCostOfTour(distanceMatrix, randomTourObj.pathOfTour);
			if (!populationSet.contains(randomTourObj)) {
				populationSet.add(randomTourObj);
				addCount++;
				Integer randomIndex = 0;
				while (true) {
					randomIndex = (int) (Math.random() * (listOfFitnessCalculatedTours.size()));
					if (randomIndex < listOfFitnessCalculatedTours.size() && randomIndex > -1
							&& !randomIndexSet.contains(randomIndex)) {
						listOfFitnessCalculatedTours.set(randomIndex, randomTourObj);
						randomIndexSet.add(randomIndex);
						break;
					}
				}
			}
		}

	}

	private ArrayList<Tour> mergeOffSpringAndOriginalPopulation(int numberOfCities, ArrayList<Tour> offSpringPopulation,
			ArrayList<Tour> listOfFitnessCalculatedTours) {
		ArrayList<Tour> mergedPopulation = new ArrayList<Tour>(listOfFitnessCalculatedTours.size());
		HashSet<ArrayList<Integer>> uniquePaths = new HashSet<ArrayList<Integer>>(listOfFitnessCalculatedTours.size());
		for (int j = 0; j < offSpringPopulation.size(); j++) {
			Tour currentOffSpring = offSpringPopulation.get(j);
			if (!uniquePaths.contains(currentOffSpring.pathOfTour)) {
				mergedPopulation.add(currentOffSpring);
				uniquePaths.add(currentOffSpring.pathOfTour);
			}
		}
		int i = 0;
		while (uniquePaths.size() != listOfFitnessCalculatedTours.size() && i < listOfFitnessCalculatedTours.size()) {
			Tour currentTour = listOfFitnessCalculatedTours.get(i);
			if (!uniquePaths.contains(currentTour.pathOfTour)) {
				mergedPopulation.add(currentTour);
				uniquePaths.add(currentTour.pathOfTour);
			}
			i++;
		}

		return mergedPopulation;
	}

	private ArrayList<Tour> crossOverSelectedPopulation(int numberOfCities,
			ArrayList<Tour> populationSelectedBasedOnFitness) {
		HashSet<ArrayList<Integer>> setOfUniqueChildren = new HashSet<ArrayList<Integer>>();
		Collections.sort(populationSelectedBasedOnFitness);
		ArrayList<Tour> offSpringPopulation = new ArrayList<Tour>();
		for (int i = 0; i < populationSelectedBasedOnFitness.size() - 1; i = i + 2) {
			ArrayList<Integer> children1 = new ArrayList<Integer>(Collections.nCopies(numberOfCities, -1));
			ArrayList<Integer> children2 = new ArrayList<Integer>(Collections.nCopies(numberOfCities, -1));
			performCrossOverOnTwoParents(numberOfCities, populationSelectedBasedOnFitness.get(i).pathOfTour,
					populationSelectedBasedOnFitness.get(i + 1).pathOfTour, children1, children2);
			if (!setOfUniqueChildren.contains(children1)) {
				setOfUniqueChildren.add(children1);
				Tour child1Tour = new Tour();
				child1Tour.pathOfTour = children1;
				offSpringPopulation.add(child1Tour);
			}
			if (!setOfUniqueChildren.contains(children2)) {
				setOfUniqueChildren.add(children2);
				Tour child2Tour = new Tour();
				child2Tour.pathOfTour = children2;
				offSpringPopulation.add(child2Tour);
			}
		}
		return offSpringPopulation;
	}

	// Using order crossover with second half as copy of parent subtour
	private void performCrossOverOnTwoParents(int numberOfCities, ArrayList<Integer> parent1,
			ArrayList<Integer> parent2, ArrayList<Integer> children1, ArrayList<Integer> children2) {
		performOrderCrossover(numberOfCities, parent1, parent2, children1);
		children1 = realignTourToCommonBase(children1);
		performOrderCrossover(numberOfCities, parent2, parent1, children2);
		children2 = realignTourToCommonBase(children2);
	}

	private void performOrderCrossover(int numberOfCities, ArrayList<Integer> parent1, ArrayList<Integer> parent2,
			ArrayList<Integer> children) {
		int midPoint = parent1.size() / 2;
		// Copy second half of parent 1
		for (int i = midPoint; i < numberOfCities; i++) {
			children.set(i, parent1.get(i));
		}
		// Copy first half in ordered manner from parent 2
		int index = 0;
		for (int i = 0; i < numberOfCities; i++) {
			int currentCity = parent2.get(i);
			if (!children.contains(currentCity)) {
				children.set(index, currentCity);
				index++;
			}
		}
	}

	private ArrayList<Tour> selectPopulationBasedOnFitnessAndRandomness(int selectedPopulationSize,
			ArrayList<Tour> listOfFitnessCalculatedTours) {
		HashSet<Tour> selectedPopulation = new HashSet<Tour>(listOfFitnessCalculatedTours.size() / 2);
		int numberOfPopulationSelectedBasedOnFitness = selectedPopulationSize / 2;
		for (int i = 0; i < numberOfPopulationSelectedBasedOnFitness; i++) {
			selectedPopulation.add(listOfFitnessCalculatedTours.get(i));
		}
		int numberOfPopulationSelectedRandomly = selectedPopulationSize - selectedPopulation.size();
		int interval = selectedPopulationSize / numberOfPopulationSelectedRandomly;
		for (int i = 0; i < numberOfPopulationSelectedRandomly;) {
			int randomTourIndex = (int) ((Math.random() * (interval)) + interval * i)
					+ numberOfPopulationSelectedBasedOnFitness;
			Tour selectedTour = listOfFitnessCalculatedTours.get(randomTourIndex);
			if (!selectedPopulation.contains(selectedTour)) {
				selectedPopulation.add(selectedTour);
				i++;
			}
		}
		return new ArrayList<Tour>(selectedPopulation);
	}

	private float injectFitnessScoresToTour(float[][] distanceMatrix, int populationSize,
			ArrayList<Tour> listOfFitnessCalculatedTours) {

		float minPathCost = Float.MAX_VALUE;
		ArrayList<Integer> maxFitTour = new ArrayList<Integer>();
		for (Tour eachTour : listOfFitnessCalculatedTours) {
			eachTour.pathCost = calculateCostOfTour(distanceMatrix, eachTour.pathOfTour);
			if (eachTour.pathCost < minPathCost) {
				minPathCost = eachTour.pathCost;
				maxFitTour = eachTour.pathOfTour;
			}
		}
		checkCostAndPrintTour(distanceMatrix, maxFitTour);
		return minPathCost;
	}

	private float determineNormalizingConstantForPopulation(float[][] distanceMatrix, int populationSize,
			ArrayList<Tour> listOfFitnessCalculatedTours) {
		int numberOfRandomCandidates = 1 + populationSize / 10;
		float sum = 0;
		int interval = populationSize / numberOfRandomCandidates;
		for (int i = 0; i < numberOfRandomCandidates; i++) {
			int randomTourIndex = (int) ((Math.random() * (interval)) + interval * i);
			ArrayList<Integer> randomlySelectedTourInInterval = listOfFitnessCalculatedTours
					.get(randomTourIndex).pathOfTour;
			sum = sum + calculateCostOfTour(distanceMatrix, randomlySelectedTourInInterval);
		}
		float meanTourCost = sum / numberOfRandomCandidates;
		return meanTourCost * 100;
	}

	private ArrayList<ArrayList<Integer>> generateInitialPopulationOfTours(int sizeOfToBeGenaratedPopulation,
			int numberOfCities, ArrayList<ArrayList<Integer>> accumulatedOverallTours) {
		HashSet<ArrayList<Integer>> initialPopulation = new HashSet<ArrayList<Integer>>();
		for (int i = 0; i < accumulatedOverallTours.size(); i++) {
			accumulatedOverallTours.set(i, realignTourToCommonBase(accumulatedOverallTours.get(i)));
		}
		int accumulatedSize = 0;
		if (accumulatedOverallTours != null) {
			initialPopulation.addAll(accumulatedOverallTours);
			accumulatedSize = initialPopulation.size();
		}
		ArrayList<Integer> randomPath = new ArrayList<Integer>();
		for (int i = 0; i < numberOfCities; i++) {
			randomPath.add(i);
		}
		while (initialPopulation.size() != sizeOfToBeGenaratedPopulation + accumulatedSize) {
			Collections.shuffle(randomPath);
			ArrayList<Integer> newRandPath = realignTourToCommonBase((ArrayList<Integer>) randomPath.clone());
			initialPopulation.add(newRandPath);
		}

		return new ArrayList<ArrayList<Integer>>(initialPopulation);
	}

	private ArrayList<Integer> greedyHeuristics(boolean isEuclideanDistance, int numberOfCities,
			float[][] euclideanDistanceMatrix, float[][] distanceMatrix) {
		ArrayList<Integer> tour = greedyHeuristicTSP(numberOfCities, distanceMatrix);
		checkCostAndPrintTour(distanceMatrix, tour);
		return tour;
	}

	private ArrayList<Integer> greedyHeuristicTSP(int numberOfCities, float[][] distanceMatrix) {
		ArrayList<Edge> listOfEdges = new ArrayList<Edge>();
		for (int i = 0; i < numberOfCities; i++) {
			for (int j = 0; j < numberOfCities; j++) {
				if (i != j) {
					Edge newEdge = new Edge();
					newEdge.cost = distanceMatrix[i][j];
					newEdge.src = i;
					newEdge.dest = j;
					listOfEdges.add(newEdge);
				}
			}
		}
		Collections.sort(listOfEdges);

		// Reuse edge class for simplicity to record prev and next links
		// src-> denotes prev, dest-> denotes next
		ArrayList<Edge> listOfNodes = new ArrayList<Edge>();
		for (int i = 0; i < numberOfCities; i++) {
			listOfNodes.add(new Edge());
		}

		int fullLinkedNodes = 0;
		for (Edge eachEdge : listOfEdges) {
			if (listOfNodes.get(eachEdge.src).dest == -1 && listOfNodes.get(eachEdge.dest).src == -1) {
				if (!findCycle(listOfNodes, eachEdge.dest, eachEdge.src)) {
					listOfNodes.get(eachEdge.src).dest = eachEdge.dest;
					listOfNodes.get(eachEdge.dest).src = eachEdge.src;

					if (listOfNodes.get(eachEdge.src).src != -1 && listOfNodes.get(eachEdge.src).dest != -1) {
						fullLinkedNodes++;
					}
					if (listOfNodes.get(eachEdge.dest).src != -1 && listOfNodes.get(eachEdge.dest).dest != -1) {
						fullLinkedNodes++;
					}
				}
			}
			if (fullLinkedNodes == numberOfCities - 2) {
				break;
			}
		}
		int nodeWithoutsrc = -1;
		for (int i = 0; i < numberOfCities && nodeWithoutsrc == -1; i++) {
			if (listOfNodes.get(i).src == -1) {
				nodeWithoutsrc = i;
			}
		}

		ArrayList<Integer> tour = new ArrayList<Integer>();
		int cur = nodeWithoutsrc;
		tour.add(cur);
		while (listOfNodes.get(cur).dest != -1) {
			cur = listOfNodes.get(cur).dest;
			tour.add(cur);
		}
		return tour;
	}

	private boolean findCycle(ArrayList<Edge> listOfNodes, int src, int dest) {
		int cur = src;
		while (listOfNodes.get(cur).dest != -1) {
			int nextNode = listOfNodes.get(cur).dest;
			if (nextNode == dest) {
				return true;
			}
			cur = nextNode;
		}
		return false;
	}

	private ArrayList<ArrayList<Integer>> savingsHeuristics(boolean isEuclideanDistance, int numberOfCities,
			float[][] euclideanDistanceMatrix, float[][] distanceMatrix) {
		ArrayList<ArrayList<Integer>> listOfTours = new ArrayList<ArrayList<Integer>>();
		int interval = 5;
		if (numberOfCities < 50) {
			interval = 1;
		}
		int numberOfSavingsHeuristicCalls = numberOfCities / interval;
		for (int i = 0; i < numberOfSavingsHeuristicCalls; i++) {
			int startingNodeNumber = (int) ((Math.random() * (interval)) + interval * i);
			ArrayList<Integer> tour = savingHeuristicTSP(startingNodeNumber, numberOfCities, distanceMatrix);
			listOfTours.add(tour);
			checkCostAndPrintTour(distanceMatrix, tour);
		}
		return listOfTours;
	}

	private ArrayList<Integer> savingHeuristicTSP(int baseNodeNumber, int numberOfCities, float[][] distanceMatrix) {
		ArrayList<Integer> tour = new ArrayList<Integer>();
		tour.add(baseNodeNumber);
		Set<Integer> unvisitedNodeSet = new HashSet<Integer>(numberOfCities);
		for (int i = 0; i < numberOfCities; i++) {
			if (i != baseNodeNumber) {
				unvisitedNodeSet.add(i);
			}
		}
		Integer[] arrayNumbers = unvisitedNodeSet.toArray(new Integer[unvisitedNodeSet.size()]);
		Random rndm = new Random();
		int rndmNumber = rndm.nextInt(unvisitedNodeSet.size());
		int currentNodeNumber = arrayNumbers[rndmNumber];
		tour.add(currentNodeNumber);
		unvisitedNodeSet.remove(currentNodeNumber);
		for (int count = 0; count < numberOfCities - 2; count++) {
			ArrayList<Integer> unvisitedCities = new ArrayList<Integer>(unvisitedNodeSet);
			Collections.shuffle(unvisitedCities);
			boolean foundLeast = false;
			for (int i : unvisitedCities) {
				if (distanceMatrix[currentNodeNumber][i] < (distanceMatrix[currentNodeNumber][baseNodeNumber]
						+ distanceMatrix[baseNodeNumber][i])) {
					currentNodeNumber = i;
					foundLeast = true;
					break;
				}
			}
			if (!foundLeast) {
				currentNodeNumber = unvisitedCities.get(0);
			}
			unvisitedNodeSet.remove(currentNodeNumber);
			tour.add(currentNodeNumber);
		}

		return tour;

	}

	private ArrayList<ArrayList<Integer>> nearestNeighbour(boolean isEucledian, int numberOfCities,
			float[][] euclideanDistanceMatrix, float[][] distanceMatrix) {
		ArrayList<ArrayList<Integer>> listOfTours = new ArrayList<ArrayList<Integer>>();
		int interval = 2;
		if (numberOfCities < 50) {
			interval = 1;
		}
		int numberOfRandomNearestNeighbourCalls = numberOfCities / interval;
		for (int i = 0; i < numberOfRandomNearestNeighbourCalls; i++) {
			int startingNodeNumber = (int) ((Math.random() * (interval)) + interval * i);
			if (isEucledian) {
				ArrayList<Integer> tour = nearestNeighbourTSP(startingNodeNumber, numberOfCities, distanceMatrix);
				listOfTours.add(tour);
				checkCostAndPrintTour(distanceMatrix, tour);
			} else {
				ArrayList<Integer> tour = nearestNeighbourTSP(startingNodeNumber, numberOfCities, distanceMatrix);
				listOfTours.add(tour);
				checkCostAndPrintTour(distanceMatrix, tour);
				tour = nearestNeighbourTSP(startingNodeNumber, numberOfCities, euclideanDistanceMatrix);
				listOfTours.add(tour);
				checkCostAndPrintTour(distanceMatrix, tour);
			}
		}
		return listOfTours;
	}

	private void checkCostAndPrintTour(float[][] distanceMatrix, ArrayList<Integer> tour) {
		float currentTourCost = calculateCostOfTour(distanceMatrix, tour);
		if (currentTourCost < leastTourCost) {
			printTour(tour);
			leastTourCost = currentTourCost;
			leastCostTour = tour;
		}
	}

	private void printTour(ArrayList<Integer> tour) {
		String tourStr = "";
		for (Integer eachCityInTour : tour) {
			tourStr = tourStr + " " + eachCityInTour;
		}
		System.out.println(tourStr.trim());
	}

	private float calculateCostOfTour(float[][] distanceMatrix, ArrayList<Integer> tour) {
		float cost = 0;
		int startingPoint = tour.get(0);
		int currentPoint = startingPoint;
		int nextPoint = startingPoint;
		for (int i = 1; i < tour.size(); i++) {
			nextPoint = tour.get(i);
			cost = cost + distanceMatrix[currentPoint][nextPoint];
			currentPoint = nextPoint;
		}
		cost = cost + distanceMatrix[currentPoint][startingPoint];
		return cost;
	}

	private ArrayList<Integer> nearestNeighbourTSP(int startingNodeNumber, int numberOfCities,
			float[][] consideredDistanceMatrix) {
		Set<Integer> visitedNodeNumber = new HashSet<Integer>();
		visitedNodeNumber.add(startingNodeNumber);
		int currentVisitedNodeNumber = startingNodeNumber;
		ArrayList<Integer> tour = new ArrayList<Integer>();
		tour.add(startingNodeNumber);
		float closestDistance = Float.MAX_VALUE;
		int nearestNeighbourCityNumber = -1;
		while (visitedNodeNumber.size() != numberOfCities) {
			closestDistance = Float.MAX_VALUE;
			for (int i = 0; i < numberOfCities; i++) {
				if (!visitedNodeNumber.contains(i)) {
					if (consideredDistanceMatrix[currentVisitedNodeNumber][i] < closestDistance) {
						closestDistance = consideredDistanceMatrix[currentVisitedNodeNumber][i];
						nearestNeighbourCityNumber = i;
					}
				}
			}
			tour.add(nearestNeighbourCityNumber);
			visitedNodeNumber.add(nearestNeighbourCityNumber);
		}
		return tour;
	}

}
