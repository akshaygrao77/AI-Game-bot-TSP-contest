package root;

import java.util.ArrayList;

public class Tour implements Comparable<Tour> {
	ArrayList<Integer> pathOfTour;
	// float fitness = -1;
	float pathCost = -1;

	public int compareTo(Tour compareTour) {
		Float pre = this.pathCost;
		Float cur = compareTour.pathCost;
		return pre.compareTo(cur);
	}
}
