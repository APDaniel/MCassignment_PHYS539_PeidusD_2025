
using MCassignment_PHYS539_PeidusD_2025.Containers;
using MCassignment_PHYS539_PeidusD_2025.Tasks;

class MyWifeHasNeverSeenMeCryingBefore
{
    public static StatisticsToExport statistics = new StatisticsToExport();
    static void Main(string[] args)
    {
        /************************************************************************************
        The script starts with a simple user's input: 1,2 or 3.
        I did not have time to design a nice interface, 
        so a desigion has been made to stick to something simple
        in order to keep it user-friendy as much as possible.
         ************************************************************************************/

        Console.WriteLine("PHYS-539 Assignment 2025. By Peidus Daniel.");
        Console.WriteLine("==========================================");
        Console.WriteLine("Choose Simulation");
        Console.WriteLine("1. Photon Kernel Simulation");
        Console.WriteLine("2. Broad paralell photon beam simulation");
        Console.WriteLine("3. Broad paralell electron beam simulation");
        Console.Write("Enter your choice (1, 2, 3): ");

        /************************************************************************************
        This part commuticates with the selected controller for each particulat simulation.
        3 simulations modeled: 
        - Monte Carlo simulations: for photon pencil beam (2,6,10 MeV); 
        for beta- field (10MeV, FS10x10);
        - Convolution-superposition: 10x10 photon field. Beam broadening effect has been omitted. 
        /************************************************************************************/

        int problem = int.Parse(Console.ReadLine());
        if (problem == 1) { RunProblemI(statistics); }
        else if (problem == 2) { RunProblemIa(statistics); }
        else if (problem == 3) { RunProblemII(statistics); }
        else Console.WriteLine("Invalid Selection");
    }
    static void RunProblemI(StatisticsToExport statistics)
    {
        Console.WriteLine("\nProblem I - Monte Carlo Simulation of Photons.");

        Problem1.KernelGenerationTask(statistics);
    }
    static void RunProblemIa(StatisticsToExport statistics)
    {
        Console.WriteLine("\nProblem Ia - Convolution-Superposition task.");

        Problem1a.InterpolationSuperpositionTask(statistics);
    }
    static void RunProblemII(StatisticsToExport statistics)
    {
        Console.WriteLine("\nProblem II - Monte Carlo Simulation of a broad Electron beam.");

        Problem2.ElectromBroadBeamTask(statistics);
    }
}