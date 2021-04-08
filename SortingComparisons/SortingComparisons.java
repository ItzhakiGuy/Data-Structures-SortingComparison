package HW4;

import java.util.Random;

import org.jfree.chart.axis.Axis;
import org.jfree.chart.needle.MiddlePinNeedle;

import Plotter.Plotter;


public class SortingComparisons {

	final static int INSERTION_VS_QUICK_LENGTH = 12;
	final static int MERGE_VS_QUICK_LENGTH = 15;
	final static int INSERTION_VS_QUICK_SORTED_LENGTH = 12;
	final static int ARBITRARY_VS_RANDOM_LENGTH = 16;
	final static int COUNTING_VS_QUICK_LENGTH = 15;
	final static double T = 600.0;

	private static void swap (double[]arr,int a, int b){
		double temp = arr[b];
		arr[b]= arr[a];
		arr[a]=temp;
	}

	/**
	 * Sorts a given array using the quick sort algorithm.
	 * At each stage the pivot is chosen to be the rightmost element of the subarray.
	 * 
	 * Should run in average complexity of O(nlog(n)), and worst case complexity of O(n^2)
	 * 
	 * @param arr - the array to be sorted
	 */
	public static void quickSortArbitraryPivot(double[] arr){
		quickSortArbitrary(arr, 0, arr.length-1);
	}
	private static void quickSortArbitrary (double[] arr, int low, int high){
		if(low<high){
			int partitionIndex = partitionArbitrary(arr, low, high);
			quickSortArbitrary(arr, low, partitionIndex-1);
			quickSortArbitrary(arr, partitionIndex+1, high);
		}
	}
	private static int partitionArbitrary (double[]arr,int low,int high){
		double pivot = arr[high];
		int i = low-1;
		for (int j = low; j < high; j++) {
			if (arr[j]<=pivot) {
				i++;
				swap(arr, i , j);
			}
		}
		swap(arr, high, i+1);
		return i+1;
	}

	/**
	 * Sorts a given array using the quick sort algorithm.
	 * At each stage the pivot is chosen in the following way:
	 * Choose a random index from the range, the element at this index is the pivot.
	 * 
	 * Should run in average complexity of O(nlog(n)), and worst case complexity of O(n^2)
	 * 
	 * @param arr - the array to be sorted
	 */
	public static void quickSortRandomPivot(double[] arr){
		quickSortRandom(arr, 0, arr.length-1);
	}
	private static void quickSortRandom (double[] arr,int low,int high){
		if(low<high){
			int p = (int)((high-low+1)*(Math.random())+low);
			swap(arr, p, high);
			int partitionIndex = partitionArbitrary(arr, low, high);
			quickSortRandom(arr, low, partitionIndex-1);
			quickSortRandom(arr, partitionIndex+1, high);
		}
	}

	/**
	 * Sorts a given array using the merge sort algorithm.
	 * 
	 * Should run in complexity O(nlog(n)) in the worst case.
	 * 
	 * @param arr - the array to be sorted
	 */
	public static void mergeSort(double[] arr)
	{
		mergeSort(arr,0,arr.length-1);
	}
	private static void mergeSort(double[] arr,int low,int high)
	{
		if (low < high)
		{
			int mid = (high+low)/2;
			mergeSort(arr, low, mid);
			mergeSort(arr, mid+1, high);
			merge(arr,low,mid,high);
		}
		else // there is only 1 item left in the sub-list
		{      
			return;                      
		}
	}
	private static void merge(double arr[], int low, int mid, int high) 
	{
        double temp1[] = new double [mid - low + 1]; 
        double temp2[] = new double [high - mid]; 
        for (int i=0; i<(mid - low + 1); i++)
        {
            temp1[i] = arr[low + i]; 
        }
        for (int j=0; j<high - mid; j++)
        {
            temp2[j] = arr[mid + 1+ j]; 
        }
        int i = 0; 
        int j = 0; 	  
        int k = low; 
        while (i < mid - low + 1 && j < high - mid) 
        { 
            if (temp1[i] <= temp2[j]) 
            { 
                arr[k] = temp1[i]; 
                i++; 
            } 
            else
            { 
                arr[k] = temp2[j]; 
                j++; 
            } 
            k++; 
        } 
        while (i < mid - low + 1) 
        { 
            arr[k] = temp1[i]; 
            i++; 
            k++; 
        } 
        while (j < high - mid) 
        { 
            arr[k] = temp2[j]; 
            j++; 
            k++; 
        } 
    } 
	 
	/**
	 * Sorts a given array, using the counting sort algorithm.
	 * You may assume that all elements in the array are between 0 and k (not including k).
	 * 
	 * Should run in complexity O(n + k) in the worst case.
	 * 
	 * @param arr
	 * @param k - an upper bound for all elements in the array.
	 */
	
	/** (stable counting sort)
	 
	public static void countingSort2(int[] arr, int k)
	{
		int[] c= new int [k];
		for (int i = 0; i < arr.length; i++)
		{
			c[arr[i]]++;
		}
		for (int i = 0; i < c.length-1; i++)
		{
			c[i+1]= (c[i])+(c[i+1]);
		}
		int [] b = new int [arr.length];
		for (int i = arr.length-1; i >= 0; i--)
		{
			b[c[arr[i]]-1] = arr[i];
			c[arr[i]]--;
		}
		for (int i = 0; i < arr.length; i++)
		{
			arr[i] = b[i];
		}
	}
	**/
	
	//unstable counting sort
	public static void countingSort(int[] arr, int k)
	{
		int[] c= new int [k];
		for (int i = 0; i < arr.length; i++)
		{
			c[arr[i]]++;
		}
		int pointer = 0;
		for (int i = 0; i < c.length; i++) {
			while (c[i]>0){
				arr[pointer] = i;
				pointer++;
				c[i]--;
			}
		}
			
	}	
	/**
	 * Sorts a given array using insertion sort.
	 * 
	 * The algorithm should run in complexity O(n^2) in the worst case.
	 * 
	 * @param arr - the array to be sorted
	 */
    public static void insertionSort(double[] arr)
    {  
    	int i = 1;
        while (i<arr.length)
        {
        	int j = i;
        	while  (j > 0 && arr[j-1]>arr[j])
        	{
				swap(arr, j-1, j);
				j--;
        	}
        	i++;
        }		
    }
    
	public static void main(String[] args)
	{		
		//insertionVsQuick();
		//mergeVsQuick();
		//insertionVsQuickOnSortedArray();
		//countingVsQuick();
		arbitraryPivotVsRandomPivot();
		
	}

	private static void countingVsQuick()
	{
		double[] quickTimes = new double[COUNTING_VS_QUICK_LENGTH];
		double[] countingTimes = new double[COUNTING_VS_QUICK_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < COUNTING_VS_QUICK_LENGTH; i++) {
			long sumQuick = 0;
			long sumCounting = 0;
			for(int k = 0; k < T; k++){
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				int[] b = new int[size];
				for (int j = 0; j < a.length; j++) {
					b[j] = r.nextInt(size);
					a[j] = b[j];
				}
				startTime = System.currentTimeMillis();
				quickSortArbitraryPivot(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				countingSort(b, size);
				endTime = System.currentTimeMillis();
				sumCounting += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			countingTimes[i] = sumCounting/T;
		}
		Plotter.plot("Counting sort on arrays with elements < n", countingTimes, "Quick sort on arrays with elements < n", quickTimes);
		
	}

	/**
	 * Compares the selection sort algorithm against quick sort on random arrays
	 */
	public static void insertionVsQuick()
	{
		double[] quickTimes = new double[INSERTION_VS_QUICK_LENGTH];
		double[] insertionTimes = new double[INSERTION_VS_QUICK_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < INSERTION_VS_QUICK_LENGTH; i++) {
			long sumQuick = 0;
			long sumInsertion = 0;
			for(int k = 0; k < T; k++){
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				quickSortArbitraryPivot(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				insertionSort(b);
				endTime = System.currentTimeMillis();
				sumInsertion += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			insertionTimes[i] = sumInsertion/T;
		}
		Plotter.plot("quick sort on random array", quickTimes, "insertion sort on random array", insertionTimes);
	}
	
	/**
	 * Compares the merge sort algorithm against quick sort on random arrays
	 */
	public static void mergeVsQuick()
	{
		double[] quickTimes = new double[MERGE_VS_QUICK_LENGTH];
		double[] mergeTimes = new double[MERGE_VS_QUICK_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < MERGE_VS_QUICK_LENGTH; i++) {
			long sumQuick = 0;
			long sumMerge = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				quickSortArbitraryPivot(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				mergeSort(b);
				endTime = System.currentTimeMillis();
				sumMerge += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			mergeTimes[i] = sumMerge/T;
		}
		Plotter.plot("quick sort on random array", quickTimes, "merge sort on random array", mergeTimes);
	}

	/**
	 * Compares the merge sort algorithm against quick sort on pre-sorted arrays
	 */
	public static void insertionVsQuickOnSortedArray()
	{
		double[] quickTimes = new double[INSERTION_VS_QUICK_SORTED_LENGTH];
		double[] insertionTimes = new double[INSERTION_VS_QUICK_SORTED_LENGTH];
		long startTime, endTime;
		for (int i = 0; i < INSERTION_VS_QUICK_SORTED_LENGTH; i++) {
			long sumQuick = 0;
			long sumInsertion = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = j;
					b[j] = j;
				}
				startTime = System.currentTimeMillis();
				quickSortArbitraryPivot(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				insertionSort(b);
				endTime = System.currentTimeMillis();
				sumInsertion  += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			insertionTimes[i] = sumInsertion/T;
		}
		Plotter.plot("quick sort on sorted array", quickTimes, "insertion sort on sorted array", insertionTimes);
	}

	/**
	 * Compares the quick sort algorithm once with a choice of an arbitrary pivot and once with a choice of a random pivot
	 */
	public static void arbitraryPivotVsRandomPivot()
	{
		double[] arbitraryTimes = new double[ARBITRARY_VS_RANDOM_LENGTH];
		double[] randomTimes = new double[ARBITRARY_VS_RANDOM_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < ARBITRARY_VS_RANDOM_LENGTH; i++) {
			long sumArbitrary = 0;
			long sumRandom = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				quickSortArbitraryPivot(a);
				endTime = System.currentTimeMillis();
				sumArbitrary += endTime - startTime;
				startTime = System.currentTimeMillis();
				quickSortRandomPivot(b);
				endTime = System.currentTimeMillis();
				sumRandom += endTime - startTime;
			}
			arbitraryTimes[i] = sumArbitrary/T;
			randomTimes[i] = sumRandom/T;
		}
		Plotter.plot("quick sort with an arbitrary pivot", arbitraryTimes, "quick sort with a random pivot", randomTimes);
	}
	

}
