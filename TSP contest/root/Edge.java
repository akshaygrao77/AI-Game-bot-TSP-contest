package root;

public class Edge implements Comparable<Edge> {
	int src = -1, dest = -1;
	float cost = -1;

	public int compareTo(Edge compareEdge) {
		Float pre = this.cost;
		Float cur = compareEdge.cost;
		return pre.compareTo(cur);
	}
}
