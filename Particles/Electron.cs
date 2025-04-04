﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace MCassignment_PHYS539_PeidusD_2025.Particles
{
    public class Electron
    {
        public Vector3 Position { get; set; }
        public Vector3 Direction { get; set; }
        public double Energy { get; set; }
        public bool IsAlive { get; set; }

        public Electron(Vector3 position, Vector3 direction, double energy, bool isAlive)
        {
            Position = position;
            Direction = direction;
            Energy = energy;
            IsAlive = isAlive;
        }
    }
}
