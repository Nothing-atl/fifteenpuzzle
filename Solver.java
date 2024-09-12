//Author : Aleena Treesa Leejoy
// 301453262
// atl17


package fifteenpuzzle;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.util.Scanner;
import java.util.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;


public class Solver {
	Node candidate;
	LinkedList<LinkedList<Node>> solutionCollection = new LinkedList<>();
	LinkedList<Node> answer = new LinkedList<>();
	double startTime;
	double endTime;
	
	
	// checks whether the program is run with at least two command-line arguments
	public static void main(String [] args) throws IOException {
		if (args.length < 2) {
			System.out.println("File names are not specified");
				System.out.println("usage: java " + MethodHandles.lookup().lookupClass().getName() + " input_file output_file");
				return;
		}
		File input = new File(args[0]); // Create a File object from the first command-line argument
		Solver compute = new Solver(input);
		String fileName = args[1]; // Write the output solution to the specified file
		List<String> output = compute.ansFormat();
		compute.writePuzzleSolution(output,fileName);
	}
	
	//The writePuzzleSolution method takes a List<String> of output solution and a String file name as
	//input parameters and writes the output to the specified file.
	public void writePuzzleSolution(List<String> outputSolution, String fileName) throws IOException {
		Path path = Paths.get(fileName);
		try {
				if (!Files.exists(path)) {
					Files.createFile(path);
			}
			Files.write(path, outputSolution); //writes the output solution to the file
			System.out.println("Solution in file: " + fileName);
		} catch (IOException e) {
			
			System.err.println("Error writing to file: " + e.getMessage());
		}
	}
	
	// define class PuzzleNode
	class PuzzleNode {
		Node predecessor; // Declare a variable of type Node called "predecessor"
		// Declare variables of type int to hold the depth and score
		int depth;
		int score;
		
		// Define a constructor that takes a Node and a score as arguments
		public PuzzleNode(Node predecessor, int score) {
			this.predecessor = predecessor;
			this.score = score;
		}
		// Define a constructor that takes a Node, a depth, and a score as arguments
		public PuzzleNode(Node predecessor, int depth, int score) {
			this.predecessor = predecessor;
			this.depth = depth;
			this.score = score;
		}
	}
	
	// Constructor method that takes in a file containing a puzzle
	public Solver(File in) {
		Node puzzle;
		Node solved;
		LinkedList<Node> temp;
		try {
			puzzle = new Node(in);
			solved = puzzle.getSolution();
			int diff = puzzle.getSize() - 5;  // Determining the difference between the puzzle size and the ideal size of 5
			startTime = System.currentTimeMillis();
			Node new_p = puzzle; // Creating new node objects to be used in the A* algorithm
			Node new_s = solved;
			if(diff > 0 && diff <= 2) { // for 6x6 and 7x7
				temp = aStarSolver(HeuristicType.manhattanandmax, new_p, new_s);
			}
			else if(diff > 2) { // for >7x7
				temp = aStarSolver(HeuristicType.euclideandis, new_p, new_s);
			}
			else {// 5x5 to less
				temp = aStarSolver(HeuristicType.manhatandis, new_p, new_s);
			}
			solutionCollection.add(temp);
		}
		catch (Exception e) {
			System.out.println(in + "not able to read" + e);
			e.printStackTrace();
		}
	}
	
	// types of heuristic used
	public enum HeuristicType {
		manhatandis,
		manhattanandmax,
		euclideandis,
	}
	
	public LinkedList<Node> aStarSolver(HeuristicType heuristicType, Node board, Node solution) {
		// Create a hash map to store puzzle nodes and their associated details
		Map<Node, PuzzleNode> nodes = new HashMap<>(); 
		// Create a priority queue to store the nodes to be checked
		Comparator<Node> scoreCompare = (a,b) -> nodes.get(a).score - nodes.get(b).score;
		PriorityQueue<Node> nodesToCheck = new PriorityQueue<>(10000, scoreCompare);
		// Initialize the starting node and add it to the hash map and priority queue
		nodes.put(board,new PuzzleNode(null,heuristicConstantCost(heuristicType,board,solution, 0)));
		nodesToCheck.add(board);
		
		
		int totalVisitedNodes = 0;
		while (!nodesToCheck.isEmpty()){
			// Remove the node with the lowest score from the priority queue
			candidate = nodesToCheck.remove();
			totalVisitedNodes++;
			// Check if the current node is the solution
			if (candidate.isSolved()){
				LinkedList<Node> solutionPath = new LinkedList<>();
				System.out.printf("number of nodes visited %d \n ",totalVisitedNodes);
				// If it is, backtrack through the nodes to create the solution path
				Node backTrace = candidate;
				while (backTrace != null){
					solutionPath.add(backTrace);
					backTrace = nodes.get(backTrace).predecessor;
				}
				return solutionPath;
			}
			// If the current node is not the solution, generate its adjacent nodes
			List<Node> adjacentNodes = AdjacentNodes(candidate);
			// Iterate over the adjacent nodes and add them to the hash map and priority queue if necessary
			adjacentNodes.forEach(node ->{
				if (!nodes.containsKey(node)){
					int newScore = heuristicConstantCost(heuristicType,node,solution, 0);
					nodes.put(node,new PuzzleNode(candidate, newScore));
					nodesToCheck.add(node);
				}
			});
		}
		return null;
	}
	
	// This method calculates the Hamming distance between two nodes.
	// It returns the number of tiles that are not in their correct position in the candidate puzzle
	// compared to their correct position in the solution puzzle.
	// but this hamming method is not used in the solution
	private int HHDistance(Node candidatePuzzle, Node solution) {
		int total = 0;
		for (getAxix currentPosition : candidatePuzzle.allNodePos()){
			int currentValue = candidatePuzzle.getValue(currentPosition);
			if (currentValue > 0) {
				getAxix targetPosition = solution.findPosition(currentValue);
				if(targetPosition != currentPosition) {
					total++;
				}
			}
		}
		return total;
	}
	
	// This method calculates the Manhattan and Max distance between two nodes.
	// It returns the sum of the Manhattan distance and the maximum of the horizontal and vertical distances
	// between each tile in the candidate puzzle and its corresponding tile in the solution puzzle.
	private int fullCoordinDis(Node candidatePuzzle, Node solved) {
		int total = 0;
		for (getAxix currentPosition : candidatePuzzle.allNodePos()) {
			int currentValue = candidatePuzzle.getValue(currentPosition);
			if (currentValue > 0) {
				getAxix targetPosition = solved.findPosition(currentValue);
				int horizontalDistance = Math.abs(targetPosition.x - currentPosition.x);
				int verticalDistance = Math.abs(targetPosition.y - currentPosition.y);
				total += Math.max(verticalDistance, horizontalDistance);
			}
		}
		total = total+ MHDistance(candidatePuzzle, solved);
		return total;
	}
	
	// This method calculates the Euclidean distance between two nodes.
	// It returns the square root of the sum of the squares of the horizontal and vertical distances
	// between each tile in the candidate puzzle and its corresponding tile in the solution puzzle.
	private int EHDistance(Node candidatePuzzle, Node solution) {
		int total = 0;
		for (getAxix currentPosition : candidatePuzzle.allNodePos()) {
			int currentValue = candidatePuzzle.getValue(currentPosition);
			if (currentValue > 0) {
				getAxix targetPosition = solution.findPosition(currentValue);
				int horizontalDistance = targetPosition.x - currentPosition.x;
				int verticalDistance = targetPosition.y - currentPosition.y;
				total +=(int)Math.sqrt(Math.pow(horizontalDistance,2)+Math.pow(verticalDistance,2));
			}
		}
		return total;
	}
	
	// This method calculates the Manhattan distance between two nodes.
	// It returns the sum of the horizontal and vertical distances between each tile in the candidate puzzle
	// and its corresponding tile in the solution puzzle.
	private int MHDistance(Node candidatePuzzle, Node solution) {
		int total = 0;
		for (getAxix currentPosition : candidatePuzzle.allNodePos()){
			int currentValue = candidatePuzzle.getValue(currentPosition);
			if (currentValue > 0) {
				getAxix targetPosition = solution.findPosition(currentValue);
				int horizontalDistance = Math.abs(targetPosition.x - currentPosition.x);
				int verticalDistance = Math.abs(targetPosition.y - currentPosition.y);
				total = total + verticalDistance + horizontalDistance;
			}
		}
		return total;
	}
	
	//The first function AdjacentNodes takes a Node object as input and returns a list of all 
	//adjacent Node objects.
	public List<Node> AdjacentNodes(Node currentBoard){
		ArrayList<Node> adjacentNodes = new ArrayList<>();
		ArrayList<getAxix> options = currentBoard.movementOptions();
		options.forEach(option-> {
			Node newPuzzle = new Node(currentBoard,option);  // create a new board by moving the blank space to that option
			adjacentNodes.add(newPuzzle);// add the new board to the adjacentNodes list
		});
		return adjacentNodes;
	}
	
	// heuristicConstantCost takes four arguments: the chosen heuristic type (HeuristicType),
	//the current board (candidatePuzzle), the goal board (solution), and an integer value (pos).
	public int heuristicConstantCost(HeuristicType heuristicChosen, Node candidatePuzzle , Node solution, int pos) {
		switch (heuristicChosen) {
			case euclideandis:
				return EHDistance(candidatePuzzle, solution);
			case manhatandis:
				return MHDistance(candidatePuzzle, solution);
			case manhattanandmax:
				return fullCoordinDis(candidatePuzzle, solution);
			default:
				throw new IllegalArgumentException("Invalid");
		}
	}
	
	
	
	//This method formats the solution found by the puzzle solver into a list of strings,
	//each representing a move made to reach the solution.
	//@return The formatted list of moves
	private ArrayList<String> ansFormat() {
		ArrayList<String> outputs = new ArrayList<>();
		Node current;
		Node next;
		getAxix newEmptyPosition;
		getAxix valueNewPosition;
		char direction = 'A';  // A default value to avoid compile-time errors
		int movedValue;
		int curSize;
		int nextSize;
		LinkedList<Node> temp = new LinkedList<>();
		
		// Extract the solution sequence from the collection of solutions found
		while(!solutionCollection.isEmpty()) {
			temp = solutionCollection.removeLast();
			while(!temp.isEmpty()) {
				answer.add(temp.removeFirst());
			}
		}
		
		// Calculate the total number of moves and the time taken to solve the puzzle
		int totalMoves = answer.size() -1;
		current = answer.removeLast();
		curSize = current.getSize();
		endTime = System.currentTimeMillis();
		double time = (endTime - startTime) / 1000.00;
		System.out.println("Time :" + time + " secs.");
		System.out.println("moves taken: " + totalMoves);
		
		// Iterate over the solution sequence to generate the list of moves
		while(!answer.isEmpty()) {
			next = answer.removeLast();
			nextSize = next.getSize();
			if(curSize != nextSize) {
				current = next;
				curSize = nextSize;
			}
			else {
				// Find the empty cell position and the moved value in the current and next nodes
				newEmptyPosition = next.findPosition(0);
				movedValue = current.getValue(newEmptyPosition);
				valueNewPosition = next.findPosition(movedValue);
				// Determine the direction of the move and add the move to the list of outputs
				if((valueNewPosition.getX() > newEmptyPosition.getX())) {
					direction = 'R';
				}
				else if((valueNewPosition.getX() < newEmptyPosition.getX())) {
					direction = 'L';
				}
				else if((valueNewPosition.getY() < newEmptyPosition.getY())) {
					direction = 'U';
				}
				else {
					direction = 'D';
				}
				String in = movedValue + " " + direction;
				outputs.add(in);
				// Update the current node and its size
				current = next;
				curSize = nextSize;
			}
		}
		return outputs;
	}
	
	
	//defines a public static class getAxis with two public fields x and y of type int
	public static class getAxix {
		public int x;
		public int y;
		
		// Define constructor that takes y and x parameters and initializes x and y fields
		public getAxix(int y, int x) {
			this.x = x;
			this.y = y;
		}
		// Define method that returns the value of x field
		public int getX() {
			return this.x;
		}
		// Define method that returns the value of y field
		public int getY() {
			return this.y;
		}
	}
	
	// creating a class Node
	public static class Node {
		private final int size;
		private int [][] board;
		private getAxix empty;
		private Node solution;
		File BoardFile;
		Scanner setup;
		public Node(File input) throws IOException {
			int y = 0;
			this.BoardFile = input;
			this.setup = new Scanner(this.BoardFile);
			String lineIn = this.setup.nextLine();
			size = Integer.parseInt(lineIn);
			board = new int [size][size];
			// Read board input from file
			while (this.setup.hasNext())
			{
				lineIn = this.setup.nextLine();
				for(int x = 0; x < size; x++) {
					try {
						// Parse integer value from substring of lineIn
						this.board[y][x] = Integer.parseInt(lineIn.substring(3 * x, 2 + 3 * x).trim());
					}
					catch (NumberFormatException e) {
						// If substring cannot be parsed, set value to 0 and store empty cell coordinates
						empty = new getAxix(y,x);
						this.board[y][x] = 0;
					}
				}
				y++;
			}
			if (BoardFile == null){
				throw new FileNotFoundException("File Not Found Exception: input file could not be found");
			}
			setup.close();
			solution = new Node(size);
		}
		
		// Constructor for creating a board with ordered tiles
		public Node(int sized) {
			size = sized;
			int counter = 1;
			board = new int [size][size];
			// Initialize the empty tile at the bottom right corner
			empty = new getAxix((size-1), (size-1));
			board[empty.getY()][empty.getX()] = 0;
			for (int y = 0; y < size; y++) {
				for (int x = 0; x < size; x++) {
					// Add numbers 1 to size^2-1 to the board in order
					board[y][x] = counter;
					counter++;
				}
			}
		}
		
		// Constructor for creating a board by moving a tile in another board
		public Node(Node currentBoard, getAxix candidate) {
			size = currentBoard.size;
			board = new int [size][size];
			// Copy the numbers from the current board to the new board
			for(int y = 0; y < size; y++) {
				for(int x = 0; x < size; x++) {
					board[y][x] = currentBoard.board[y][x];
				}
			}
			// Copy the location of the empty tile
			empty = currentBoard.getEmpty();
			if(isValid(candidate)) {
				tileMove(candidate);
			}
			solution = new Node(size);
		}
		
		// Getter for the solution board
		public Node getSolution() {
			return this.solution;
		}
		
		// Getter for the location of the empty tile
		private getAxix getEmpty() {
			return empty;
		}
		
		// Getter for the value at a certain location
		public int getValue(getAxix location) {
			return board[location.getY()][location.getX()];
		}
		
		// Getter for the size of the board
		public int getSize() {
			return size;
		}
		
		// Checks if a candidate move is valid
		private boolean isValid(getAxix candidate) {
			if((candidate.getY() < 0||candidate.getY() >= size) || (candidate.getX() < 0||candidate.getX() >= size)) {
				return false;
			}
			if((Math.abs(empty.getY() - candidate.getY()) + Math.abs(empty.getX() - candidate.getX())) != 1) {
				return false;
			}
			return true;
		}
		
		// Checks if a candidate move is valid
		public void tileMove(getAxix location) throws IllegalArgumentException{
			if(!isValid(location)) {
				throw new IllegalArgumentException("Illegal Argument Exception: The desired movement is not possible");
			}
			board[empty.getY()][empty.getX()] = board[location.getY()][location.getX()];
			empty = location;
			board[location.getY()][location.getX()] = 0;
		}
		
		// Finds all possible movement options for the empty space
		public ArrayList<getAxix> movementOptions() {
			ArrayList<getAxix> options = new ArrayList<>();
			getAxix candidate;
			int x;
			int y;
			for(int shift_y = -1; shift_y <= 1 ; shift_y++) {
				for(int shift_x = -1; shift_x <= 1; shift_x++) {
					x = empty.getX() + shift_x;
					y = empty.getY() + shift_y;
					candidate = new getAxix(y,x);
					if(isValid(candidate)) {
						options.add(candidate);
					}
				}
			}
			return options;
		}
		
		// Finds positions of all nodes on the board
		public List<getAxix> allNodePos() {
			ArrayList<getAxix> movements = new ArrayList<>();
			for(int y = 0; y < size; y++) {
				for(int x = 0; x < size; x++) {
					movements.add(new getAxix(y,x));
				}
			}
			return movements;
		}
		
		// Finds the position of a specific item on the board
		public getAxix findPosition(int item) {
			getAxix output;
			for(int y = 0; y < size; y++) {
				for(int x = 0; x < size; x++) {
					
					if(board[y][x] == item) {
						output = new getAxix(y,x);
						return output;
					}
				}
			}
		   return null;
		}
		
		
		 // Check if the current board configuration is solved.
		 //@return boolean value indicating if the board is solved or not.
		public boolean isSolved() {
			for(int y = 0; y < size; y++) {
				for(int x = 0; x < size; x++) {
					
					if((board[y][x] > 0) && (board[y][x] != solution.board[y][x])) {
						return false;
					}
				}
			}
			return true;
		}
		
		//Override the equals method to compare two objects.
		//o Object to be compared with this Node instance.
		@Override
		public boolean equals(Object o) {
			if (this == o) {
				return true;
			}
			if (o == null || getClass() != o.getClass()) {
				return false;
			}
			Node other = (Node) o;
			for (int y = 0; y < size; y++) {
				for (int x = 0; x < size; x++) {
					if (this.board[y][x] != other.board[y][x]) {
						return false;
					}
				}
			}
			return true;
		}
		
		//compute a hash code for the Node object.
		//return an integer value representing the hash code of this Node instance.
		@Override
		public int hashCode() {
			int result = 17;
			for (int y = 0; y < size; y++) {
				for (int x = 0; x < size; x++) {
					result = 31 * result + board[y][x];
				}
			}
			return result;
		}
	}
}