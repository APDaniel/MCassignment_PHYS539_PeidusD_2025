
using MCassignment_PHYS539_PeidusD_2025.Tasks;

class MyWifeHasNeverSeenMeCryingBefore
{
    static void Main(string[] args)
    {
        Console.WriteLine("PHYS-539 Assignment 2025. By Peidus Daniel.");
        Console.WriteLine("==========================================");
        Console.WriteLine("Choose Problem");
        Console.WriteLine("1. Problem I  (Photon Simulation)");
        Console.WriteLine("2. Problem II (Electron Simulation)");
        Console.Write("Enter your choice (1 or 2): ");

        int problem = int.Parse(Console.ReadLine());
        if (problem == 1) { RunProblemI(); }
        else if (problem == 2) { RunProblemII(); }
        else Console.WriteLine("Invalid Selection");
    }
    static void RunProblemI()
    {
        Console.WriteLine("\nProblem I - Monte Carlo Simulation of Photons.");
        Console.WriteLine("Select task:");
        Console.WriteLine("1. Generate primary dose deposition kernel.");
        Console.WriteLine("2. Plot kernel isolines.");
        Console.WriteLine("3. Plot Compton angle distribution.");
        Console.WriteLine("4. Compute total primary energy deposited.");
        Console.WriteLine("5. Convolution for dose in water phantom.");
        Console.WriteLine("6. Plot PDD and Profiles.");
        Console.Write("Enter task number (1-6): ");
        int task = int.Parse(Console.ReadLine());
        
        switch (task)
        {
            case 1:
                Problem1Task1.KernelGenerationTask();
                break;
            default:
                Console.WriteLine("Invalid selection...");
                break;
        }
    }
    static void RunProblemII()
    {

    }
}