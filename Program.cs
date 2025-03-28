﻿
using MCassignment_PHYS539_PeidusD_2025.Containers;
using MCassignment_PHYS539_PeidusD_2025.Tasks;

class MyWifeHasNeverSeenMeCryingBefore
{
    public static StatisticsToExport statistics = new StatisticsToExport();
    static void Main(string[] args)
    {
        Console.WriteLine("PHYS-539 Assignment 2025. By Peidus Daniel.");
        Console.WriteLine("==========================================");
        Console.WriteLine("Choose Simulation");
        Console.WriteLine("1. Photon Kernel Simulation");
        Console.WriteLine("2. Broad paralell photon beam simulation");
        Console.WriteLine("3. Broad paralell electron beam simulation");
        Console.Write("Enter your choice (1 or 2): ");

        int problem = int.Parse(Console.ReadLine());
        if (problem == 1) { RunProblemI(statistics); }
        else if (problem == 2) { RunProblemII(); }
        else if (problem == 3) { RunProblemIII(); }
        else Console.WriteLine("Invalid Selection");
    }
    static void RunProblemI(StatisticsToExport statistics)
    {
        Console.WriteLine("\nProblem I - Monte Carlo Simulation of Photons.");

        Problem1.KernelGenerationTask(statistics);
    }
    static void RunProblemII()
    {

    }
    static void RunProblemIII()
    {

    }
}