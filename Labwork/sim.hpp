#ifndef _SIMCLASS_
#define _SIMCLASS_
/*
 * X = done 
 * * = todo
 * E = exercise
 * OE = optional exercise
 *
 *  ===== base class (common stuff for MC and MD) =====  [ DD/MM/2025 ] 
 *
 *   [X]  calcenergyi()
 *   [X]  calctotenergy()
 *   [X]  pbc() [ to apply pbc ]
 *   [X]  save_mgl_snapshot() save configurations in mgl format
 *   [X]  prepare_initial_conf() 
 *   [X]  init_rng()
 *   [X]  run() method (main loop) 
 *
 *  ===== MC class =====
 *
 *  NTV ensemble [ 25/11/2025 ]
 *
 *   [X]  alpha()
 *   [X]  acc()
 *   [X]  move_NTV()
 *
 *   [X]  calc_acceptance and adjust()
 *   [X]  init_measures() truncate measure files
 *   [X]  save_measures() append measures
 *
 *   [OE]  calculate tail corrections to energy and pressure
 *
 *  optional exercises:
 *
 *   // in "sim" base class
 *   [OE]  save_snapshot() save configurations 
 *   [OE]  save_restart()
 *   
 *   NPT ensemble [ 09/12/2025 ]
 *
 *   // in "mcsim" derived class
 *
 *   [EX] init_measures() and save_measures(): add density saving
 *
 *   [EX] modify run() method to attempt box moves too
 *   
 *   [ 09/12/2025 ] EXERCISES FOR NPT ensemble:
 *
 *   [EX] restore_all_pars (restore positions of all particles)
 *
 *   [EX] store_all_pars (store position of all particles)
 *
 *   [EX] alpha_box() i.e. trial box move
 *   
 *   [EX] acc_box() i.e. acceptance of trial box move
 *
 *** MD class [ 09/01/2026 ]
 *
 *   [E]  init_measures() truncate measure files
 *
 *   [E]  calcK() calculate the total kinetic energy
 *   
 *   [E]  save_measures() save measures
 *   
 *   [E]  gauss() (gauss stochastic variable using Box-Muller)
 *   
 *   [E]  calc_forces()  <---
 *   
 *   [E]  run() (implement velocity verlet as a symmetric factorization of propagator) <---
 *   
 *   [E]  prepare_initial_conf() (initialize velocities) <---
 *   
 *
 */
#include <vector>
#include <string>
#include <fstream>
#include "./params.hpp"
#include "./particle.hpp"
#include "./randnumgen.hpp"
#include <iomanip> // for setprecision()

template<typename ntype, typename particle_type>
class sim
{
  using simp = simpars<ntype>;
protected:
  simp pars;  // parametri per composizione
  std::vector<particle_type> parts;
  // if opt=1 calculate energies only for i < j
  ntype calcenergyi(int i, int opt=0) 
    {
      int j;
      ntype enei=0.0;
      for (j=0; j < pars.Np; j++) // pars.Np è il numero totale di particelle
        {
          if (opt==1 && i >= j)
            continue;
          if (i==j)
            continue;
          // la classe particelle deve implementare un metodo vij per il calcolo dell'energia d'interazione
          enei += parts[i].vij(parts[j], pars.L); 
          // pars.L è un vettore con i lati del box 
        }
      return enei;
    }
  
  ntype totenergy(void)
    {
      ntype ene=0.0;

      for (auto i=0; i < static_cast<int>(parts.size()); i++)
        {
          ene+=calcenergyi(i, 1);
        }
      return ene;
    }

  void pbc(int i) 
    {
      auto Dr = parts[i].r;
      Dr = mulcw(pars.L,rint(divcw(Dr,pars.L))); // L*rint(Dr/L)
      parts[i].r = parts[i].r - Dr;
    }

  void save_mgl_snapshot(long int t) 
    {
       std::fstream f;
       std::string s;

       s = "cnf-" + std::to_string(t) + ".mgl";
       f.open(s, std::ios::out|std::ios::trunc);
       for (int i=0; i < pars.Np; i++)
         {
           f << parts[i].r(0) << " " << parts[i].r(1) << " " <<
             parts[i].r(2) << " @ " << pars.sigma*0.5 << "\n";
         }
       f.close();
    }
 
public:
  void prepare_initial_conf(void) 
    {
      // SC 
      int ix, iy, iz;
      int cc=0;
      pars.Np = pars.nx*pars.ny*pars.nz;
      parts.resize(pars.Np);
      ntype vcell = pow(pars.sigma,3.0);
      ntype rhomax = 1.0/vcell;
      ntype sf;
      sf = cbrt(rhomax/pars.rho);
      pars.L ={ntype(pars.nx), ntype(pars.ny), ntype(pars.nz)};
      std::cout << "Box size is " << pars.L << "\n";
      std::cout << "Initial density is " << pars.Np/(pars.L(0)*pars.L(1)*pars.L(2)) << "\n";
      ntype clen = sf*pars.sigma;
      pars.L *= clen;
      for (ix = 0; ix < pars.nx; ix++)
        for (iy = 0; iy < pars.ny; iy++)
          for (iz = 0; iz < pars.nz; iz++)
            {
              parts[cc].r = {ix*clen, iy*clen, iz*clen};
              parts[cc].r -= pars.L*0.5;
              parts[cc].set_sigma(pars.sigma);
              parts[cc].set_epsilon(pars.epsilon);
              parts[cc].set_rcut(pars.rc);
              cc++;
            }
      // ...or BCC or FCC lattice
 
    }

  void set_sim_type(int type)
    {
      pars.simtype = type;
    }
 
  void init_rng(int n) 
    {
      if (n < 0)
        rng.rseed();
      else
        rng.seed(n);
    }

  void run(void) 
    {
      // intentionally void
    }; 
};

template<typename ntype, typename particle_type>
class mcsim: public sim<ntype, particle_type>
{
  // for calc_acceptance_and_adjust: total trial moves and accepted ones
  // for calculating acceptance rates.
  using bc=sim<ntype, particle_type>;
  using bc::calcenergyi, bc::parts, bc::pars, 
        bc::totenergy, bc::save_mgl_snapshot, bc::pbc;
  
  // counters used to calculate acceptance rates
  long int tot_tra, tot_vol, tra_rej, vol_rej;

  void alpha(int i)
    {
      // displacement pars.deltra
      // move particle i by using as max displacemente paramater pars.deltra
      // and apply PBC 
      pvector<ntype,3> delr;
      delr ={pars.deltra*2.0*(rng.ranf()-0.5),  pars.deltra*2.0*(rng.ranf()-0.5),
           pars.deltra*2.0*(rng.ranf()-0.5)};
      parts[i].tra_move(delr);
      pbc(i);
    }

  void acc(int i, ntype eno)
    {
      // ACCEPT OR REJECT TRIAL MOVE ACCORDING TO METROPOLIS CRITERIUM
      // We assume kB=1 (Boltzmann constants) so that beta=1/T. 
      // Calculate new energy and accept move with probability
      // min {1, exp(-\beta deltaU}, i.e.
      // if deltaU < 0 accept the move
      // otherwise generate a random number xi in [0,1] 
      // and if xi < exp(-beta*deltaU) accept the move otherwise reject it 
      ntype enn = calcenergyi(i);
      ntype delu= enn-eno;
      ntype xi = rng.ranf();
      if (delu > 0.0 && xi >= exp(-delu/pars.T))
        {
          // reject move
          tra_rej++;
          parts[i].restore();
        }
    }
 
  void move_NTV(int i)
    {
      // 1) calculate energy of particle i
      ntype eno;
      eno = calcenergyi(i);
      // 2) store position of particle i
      parts[i].store();
      // 3) trial move 
      alpha(i);
      // 4) acceptance 
      acc(i, eno);
    }

  void move_box(void)// <------------------------------------------------
    {
      ntype DG, fact; // DG is \Delta G as discussed during lectures (see pdf of lectures)
      store_all_pars(); // store all particle positions before attempting a box move
      alpha_box(DG, fact);
      acc_box(DG, fact);
    }

  void alpha_box(ntype& DG, ntype& fact)
    {
      // trial box move
      //
      // 1) choose randomly a new volume
      //
      // 2) scale all particle positions and box size by factor "fact"
      // 
      // 3) calculate new total interaction energy
      //
      // 4) calculate \DeltaG (see pdf of lectures)

      // EXERCISE SOL code here
      ntype beta, delV, Vo, Vn, enn, eno;
      delV = pars.delv*(rng.ranf()*2.0-1.0);
      eno = totenergy();
      Vo = pars.L(0)*pars.L(1)*pars.L(2);
      Vn = Vo + delV;
      fact = cbrt(Vn/Vo); // (Vn/Vo)^(1/3)
      for (auto i=0; i < pars.Np; i++)
        {
          parts[i].r *= fact;
        }
      pars.L *= fact;
      enn = totenergy();
      beta = 1.0/pars.T;
      DG = pars.P*delV + (enn-eno) - (1.0/beta)*pars.Np*log(Vn/Vo); 
    }

  void acc_box(ntype DG, ntype fact)
    {
      // accept or reject box move
      //
      // 1) calculate e^{-\Delta G} (see pdf of lectures)
      // 2) generate a random numner \xi and check wheter to accept the box trial move      
      // 3) if move is rejected:
      //    i)   restore all particle positions thourgh method restore_all_pars()
      //    ii)  restore box size
      //    iii) update counter vol_rej of rejected box moves for calculating acceptance ratio 

      // EXERCISE SOL: code here
      ntype beta, expDG, xi;
      beta=1.0/pars.T;
      expDG = exp(-beta*DG);
      xi = rng.ranf();
      if (xi >= expDG)
        {
          // reject 
          restore_all_pars();
          pars.L /= fact;
          vol_rej++;
        }
    }
 

 void restore_all_pars()
   {
     // restore all particle positions 
     // EXERCISE SOL: code here
     for (auto i=0; i < pars.Np; i++)
        {
          parts[i].restore();
        }
   }

 void store_all_pars()
   {
     // store all particle position 
     // EXERCISE SOL: code to store all parts pos
     for (auto i=0; i < pars.Np; i++)
        {
          parts[i].store();
        }
   }

  void calc_acceptance_and_adjust(void)
    {
      // CALC ACCEPTANCE RATES AND 
      // ADJUST pars.deltra 
      ntype r_tra, r_vol;
      if (tot_tra > 0)
        {
          // estimate acceptance rate
          r_tra = ((double)(tot_tra - tra_rej))/tot_tra;
          std::cout << "rate tra: " << r_tra << " deltra=" << pars.deltra << "\n";
          if (r_tra > 0.5)
            {
              pars.deltra *= 1.1;
            }
          else
            {
              pars.deltra /= 1.1;
            }
          tot_tra = tra_rej = 0;
        }

      // adjust maximum volume "displacement" in alpha_box
      // so that acceptance rate of box move is around 0.5
      if (tot_vol > 0) // <------------------------------------------------
        {
          r_vol = ((ntype)(tot_vol - vol_rej)) / tot_vol;
          std::cout << "rate vol: " << r_vol << " delv=" << pars.delv << "\n";
          if (r_vol > 0.5) 
            {
              pars.delv *= 1.1;
           }
          else
            {
              pars.delv /= 1.1;
            }
          tot_vol=vol_rej=0.0;
        }
    }

  void init_measures(void)
    {
      // open files in writing mode to truncate them to 0
      std::fstream f;
      f.open("energy.dat", std::ios::out|std::ios::trunc);
      f.close();
      if (pars.simtype==1) // simtype == 1 means an NPT simluations
        {
          // EXERCISE SOL: code to save density
          f.open("density.dat", std::ios::out|std::ios::trunc);
          f.close();
        }
      // init your own measures here
    }

 void save_measures(long int t)
    {
      std::fstream f;
      f.open("energy.dat", std::ios::out|std::ios::app);
      // save potential energy per particle
      f << t << " " << totenergy()/pars.Np << "\n";
      f.close();
      if (pars.simtype==1) // 1 means NPT simulation, save density in this case
        {
          // EXERCISE SOL: code to save density
          f.open("density.dat", std::ios::out|std::ios::app);
          f << t << " " << pars.Np/(pars.L(0)*pars.L(1)*pars.L(2)) << "\n";
          f.close();
        }
      // you can add your own measures here
    }

 public:
  void run(void) 
    {
      // loop over MC steps
      int i, t, ip;
      tot_tra = tra_rej = tra_rej = vol_rej = 0;
      init_measures();
      if (pars.simtype==0)
        std::cout << "MC-NTV\n";
      else
        std::cout << "MC-NPT\n";
      for (t = 0; t < pars.totsteps; t++)
        {
          // ATTEMPT TO MOVE ALL PARTICLES:
          // Np move attempts = 1 MC step
          for (i=0; i < pars.Np; i++)
            {
              // EXERCISE SOL: choose between box move and particle move
              if (pars.simtype==1 && rng.ranf() < 1.0/(pars.Np+1))
                {
                  move_box();// <------------------------------------------------
                  tot_vol++; 
                }
              else 
                {
                  ip = rng.ranf()*pars.Np;
                  move_NTV(ip);
                  tot_tra++;
                }
            }
        
          if (t > 0 && pars.savemeasure > 0 && t % pars.savemeasure == 0)
            {
              save_measures(t);
            }

          if (t > 0 && t % pars.outstps == 0) 
            {
              std::cout << "Step #" << t << "\n";
              // NOTE: to compare with Johnson paper data the internal energy must be calculated (i.e. by including
              // the kinetic contribution)
              std::cout << "total energy per particle is " << totenergy()/pars.Np << "\n";
            }

          if (t > 0 && pars.save_mgl_snapshot > 0 && 
              t % pars.save_mgl_snapshot == 0)
            {
              save_mgl_snapshot(t);
            }

          if (t > 0 && pars.adjstps > 0 && t % pars.adjstps == 0 && 
              t < pars.maxadjstps)
            {
              calc_acceptance_and_adjust();
            }
        }
    } 
};

template<typename ntype, typename particle_type>
class mdsim: public sim<ntype, particle_type>
{
  using bc = sim<ntype, particle_type>;
  using bc::pars, bc::parts, bc::pbc;
  ntype Us, U;
 
  void init_measures(void)
    {
      // EXERCISE write your own code here
    }

  ntype calcK(void)
    {
      // EXERCISE write your own code here
      return 0;
    }
 
  void save_measures(long int t)
    {
      // EXERCISE write your own code here
    }
  // calculate force acting on each particle and 
  // at the same time the total potential energy U
  // and the total shifted potential energy
  void calc_forces(ntype& Us, ntype& U)
    {
      // EXERCISE write your own code here
    }

  // Generate gaussian number by Box Muller algorithm  
  ntype gauss(void)
    {
      // gaussian of variance 1 and 0 mean
      double x1, x2;
      do
        {
          x1 = rng.ranf();
          x2 = rng.ranf();
        }
      while (x1==0 || x2==0);
      return cos(2*M_PI*x2)*sqrt(-2.0*log(x1));
    }

public:
  void prepare_initial_conf(void)
    {
      bc::prepare_initial_conf(); // place particles on a SC lattice

      // EXERCISE write your own code here
      // and now set velocities according to Maxwell-Boltzmann distribution
      //
      // 1) use gauss() method (Box-Muller algorithm) to initilize particle's velocity from Maxwell-Boltzmann distribution,
      // [ hint: standard deviation of particle velocity is sqrt(pars.T / parts[i].m) ]
      //
      // 2) set also vcut (i.e. v_{LJ}(r_{cut}) by calling parts[i].set_vcut() method
      //
      // 3) set total momentum to 0
      //
    }
  void run(void)
    {
      // s is the MD step
      calc_forces(Us, U);
      std::cout << "Initial Excess energy=" << Us << "\n";
      std::cout << "Initial Kinetic energy=" << calcK() << "\n";

      std::fstream fo;
      init_measures();
      for (long int t=0; t < pars.totsteps; t++)
        { 
          //////////////////////////////////////////////
          // IMPLEMENT A VELOCITY-VERLET
          //
          // exp(iL*dt)=exp(iLp*dt/2)*exp(iLq*dt)*exp(iLp*dt/2)
          //
          // WRITE THIS
          //
          // 1) Apply exp(iLp*dt/2) to all particles
          //
          // 2) Apply exp(iLq*dt) to all particles and PBC
          //
          // 3) calculate forces
          //
          // 4) apply exp(iLp*dt/2) to all particles here
          // 
          ////////////////////////////////////////////////
          // EXERCISE write your own code here
          ////////////////////////////////////////////////

          if (pars.outstps > 0 && t % pars.outstps== 0)
            {
              ntype K=calcK();
              std::cout << std::setprecision(15) << "[step #" << t << "] \nTotal energy=" << K+Us << "\n";
              std::cout << std::setprecision(8) << "Excess energy per particle=" << U/pars.Np << "\n";
              std::cout << std::setprecision(8) << "Kinetic energy=" << K << " T=" 
                << 2*K/(3*(pars.Np-3)) << "\n";
            }
          if (t > pars.eqstps && pars.savemeasure > 0 && (t % pars.savemeasure == 0))
            {
              save_measures(t);
            }

          if (t > 0 && pars.save_mgl_snapshot > 0 && 
              t % pars.save_mgl_snapshot == 0)
            {
              bc::save_mgl_snapshot(t);
            }
        }
    }
};
#endif
