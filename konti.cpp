
   /* -----------------------------------------------------------------------

      Programm:   konti
                  
      Sprache:    ANSI C++
                  
      Autor:      Martin Petry
                  
      E-Mail:     martin.petry@fh-bielefeld.de
                  
      Repository: https://github.com/cfd-bielefeld/konti-1d-implizit
                  
      Lizenz:     GNU GPL v3.0
                  
      Erstellt:   09.04.2021
                  
      Geändert:   19.04.2021

      -----------------------------------------------------------------------

      Dieses Programm löst die eindimensionale Kontinuitätsgleichung

                           rho_t + (rho u)_x = 0

      für ein gegebenes stationäres Geschwindigkeitsfeld u(x) auf einem
      Zell-zentrierten strukturierten Gitter.

      Die Berechnung der Massenströme erfolgt durch Verwendung analytischer
      Formeln zur Vermeidung von if-Abfragen.

      Da das Geschwindigkeitsfeld vorgegeben ist, wird keine Rückwirkung
      der Dichte auf die Geschwindigkeit gemäß Navier-Stokes Gleichung
      berücksichtigt. Wir können die Navier-Stokes Gleichung retten,
      indem wir geeignete Zwangskräfte einführen.

      -----------------------------------------------------------------------

      Modelleigenschaften:

           1. Eindimensional

           2. Implizit instationär

           3. Finite Volumen Methode

           4. Zell-zentriertes Gitter

      -----------------------------------------------------------------------

      Folgende Randbedingungen sind implementiert:

           1. Fester Rand
           2. Dirichlet-Bedingung
           3. Neumann Randbedingung
           4. Periodische Randbedingung
           5. Dynamische Randbedingung
           6. Strömungsausgang
      
      ----------------------------------------------------------------------- */

   /* Header Dateien einbinden */

     #include <iostream>
     #include <fstream>
     #include <iomanip>
     #include <string>
     #include <stdlib.h>
     #include <math.h>

   /* ----------------------------------------------------------------------- */

   /* Namensraum festlegen */

      using namespace std;

   /* ----------------------------------------------------------------------- */

   /* Versionsnummer definieren */

      string CppVersion = "19.04.2021";

   /* ----------------------------------------------------------------------- */

   /* Funktionsprototypen deklarieren */

      void StartmeldungAusgeben     (void);
      void ParameterEinlesen        (void);
      void SimulationVorbereiten    (void);
      void ParameterAusgeben        (void);
      void SimulationDurchfuehren   (void);
      void ErgebnisseSpeichern      (void);
      void ErgebnisseDarstellen     (void);
      void ProgrammBeenden          (void);
                                    
      void MasseBerechnen           (void);
      void EkinBerechnen            (void);
      void ImpulsBerechnen          (void);
                                    
      void (*HilfsfelderBerechnen ) (void);
      void (*ImpliziterZeitschritt) (void);


      void HilfsfelderBerechnenPBC  (void);
      void HilfsfelderBerechnenSTD  (void);
                                    
      void ImpliziterZeitschrittPBC (void);
      void ImpliziterZeitschrittSTD (void);
                                    
      void SetGradient              (double*,int,double);
      void GetGradient              (double*,int,double&);
                                    
      void RandbedingungAusgeben    (int,char*);
      void FunktionsNameAusgeben    (int,char*);
                                    
      double (*U0)                  (double);
      double (*R0)                  (double);
      double (*W0)                  (double);
      double (*O0)                  (double);
                                    
      double Konstante              (double);
      double Linear                 (double);
      double Parabel                (double);
      double Rechteck               (double);
      double Dreieck                (double);
      double Saegezahn              (double);
      double Linearrampe            (double);
      double Cosinusrampe           (double);
      double Sinus                  (double);
      double Cosinus                (double);
      double Exponential            (double);
      double Gauss                  (double);
      double Dirac                  (double);
      double Heaviside              (double);
      double CosPeak                (double);
      double UserDefined01          (double);
      double UserDefined02          (double);
                                    
      double VolumenIntegration     (double*);

   /* ------------------------------------------------------------------ */

   /* Globale Variablen deklarieren */

      double xa;     // Linke Intervallgrenze
      double xe;     // Rechte Intervallgrenze
      double x0;     // Geisterpunkt P0
      double t0;     // Globaler Anfangszeitpunkt der Simulation
      double ta;     // Lokaler Anfangszeitpunkt der Simulation
      double te;     // Endzeitpunkt der Simulation
      double dt;     // Länge des Zeitschritts
      double Pi;     // Kreiszahl Pi

      double dxrhoW; // x-Ableitung von rho auf Westrand
      double dxrhoO; // x-Ableitung von rho auf Ostrand
                    
      double u1;     // Lageparameter x-Richtung
      double u2;     // Lageparameter u-Richtung
      double u3;     // Formparameter x-Richtung
      double u4;     // Formparameter u-Richtung
      double u5;     // Summand bei xa
      double u6;     // Summand bei xe
                     
      double r1;     // Lageparameter x-Richtung
      double r2;     // Lageparameter r-Richtung
      double r3;     // Formparameter x-Richtung
      double r4;     // Formparameter r-Richtung
      double r5;     // Summand bei xa
      double r6;     // Summand bei xe
                     
      double w1;     // Lageparameter t-Richtung
      double w2;     // Lageparameter r-Richtung
      double w3;     // Formparameter t-Richtung
      double w4;     // Formparameter r-Richtung
                     
      double o1;     // Lageparameter t-Richtung
      double o2;     // Lageparameter r-Richtung
      double o3;     // Formparameter t-Richtung
      double o4;     // Formparameter r-Richtung
                    
      double t ;     // Aktuelle Zeit
                     
      int imin ;     // Index für erste Zelle
      int imax ;     // Anzahl der räumlichen Teilintervalle
      int mesh ;     // Art des Gitters ( intern / extern )
                     
      int IMAX ;     // Maximale Anzahl der Gauß-Seidel Iterationen
      int NMAX ;     // Zähler für: IMAX erreicht!
                     
      int RBW  ;     // Randbedingung am Westrand
      int RBO  ;     // Randbedingung am Ostrand
                     
      int wF   ;     // Randfunktion West für rho
      int oF   ;     // Randfunktion Ost  für rho
      int uF   ;     // Anfangsfunktion für u
      int rF   ;     // Anfangsfunktion für rho
      int ED   ;     // Ergebnisse darstellen
      int AZ   ;     // Anfangszustand
                     
      double delta;  // Absolute Genauigkeit für den Defekt

      double *x ;    // Hilfsfeld
      double *dx;    // Hilfsfeld

      double *u;      // Geschwindigkeitsfeld
      double *rho;    // Dichtefeld
                      
      double *f;      // Hilfsfeld

      double *aW;     // Nebendiagonale der Koeffizientenmatrix
      double *aP;     // Hauptdiagonale der Koeffizientenmatrix
      double *aE;     // Nebendiagonale der Koeffizientenmatrix

      unsigned long int nmax ; // Anzahl der zeitlichen Teilintervalle
      unsigned long int N    ; // Jeder N-te Zeitschritt wird gespeichert

      double M;       // Aktuelle Masse
      double Ek;      // Aktuelle kinetische Energie
      double px;      // Aktueller Impuls

      ofstream Dout;        // Ausgabeobjekt für Residuum
                            
      enum                  // Enum-Konstanten für RB
     {                      
      WallBoundary      ,   // = 0
      DirichletBoundary ,   // = 1
      NeumannBoundary   ,   // = 2
      PeriodicBoundary  ,   // = 3
      DynamicBoundary   ,   // = 4
      OutletBoundary        // = 5
     };                     
                                                              
      double (*Funktion[17]) (double) =  // Feld von Zeigern auf Funktionen
     { 
      Konstante,       //  =  0
      Linear,          //  =  1
      Parabel,         //  =  2
      Rechteck,        //  =  3
      Dreieck,         //  =  4
      Saegezahn,       //  =  5
      Linearrampe,     //  =  6
      Cosinusrampe,    //  =  7
      Sinus,           //  =  8
      Cosinus,         //  =  9
      Exponential,     //  = 10
      Gauss,           //  = 11
      Dirac,           //  = 12
      Heaviside,       //  = 13
      CosPeak,         //  = 14
      UserDefined01,   //  = 15
      UserDefined02    //  = 16
     };

   /* -----------------------------------------------------------------------
      Anfang von main
      ----------------------------------------------------------------------- */

      int main(void)
     {

      StartmeldungAusgeben();
      ParameterEinlesen();
      SimulationVorbereiten();
      ParameterAusgeben();
      SimulationDurchfuehren();
      ErgebnisseSpeichern();
      ErgebnisseDarstellen();
      ProgrammBeenden();

      return (0);

     }

   /* -----------------------------------------------------------------------
      Ende von main
      ----------------------------------------------------------------------- */



   /* -----------------------------------------------------------------------
      Anfang von StartmeldungAusgeben
      ----------------------------------------------------------------------- */

      void StartmeldungAusgeben(void)
     {

      system("clear");

      cout << "\n                           * * * Kontinuitätsgleichung 1D * * *\n"
           << "\n                                   Version: " << CppVersion << "\n\n";

     }

   /* -----------------------------------------------------------------------
      Ende von StartmeldungAusgeben
      ----------------------------------------------------------------------- */



   /* -----------------------------------------------------------------------
      Anfang von ParameterEinlesen
      ----------------------------------------------------------------------- */

      void ParameterEinlesen(void)
     {

      ifstream fin;       // Eingabeobjekt zum Lesen aus einer Datei

      string InpVersion;  // Versionsnummer der Eingabedatei

      fin.open("input.dat");

                      fin.ignore(80,'\n');
                      fin.ignore(80,'\n');
                      fin.ignore(80,'\n');
                      fin.ignore(26,'\n'); fin >> InpVersion; fin.ignore(80,'\n');
                      fin.ignore(80,'\n');
      fin >> xa;      fin.ignore(80,'\n');
      fin >> xe;      fin.ignore(80,'\n');
                      fin.ignore(80,'\n');
      fin >> ta;      fin.ignore(80,'\n');
      fin >> dt;      fin.ignore(80,'\n');
                      fin.ignore(80,'\n');
      fin >> rF;      fin.ignore(80,'\n');
                      fin.ignore(80,'\n');
      fin >> r1;      fin.ignore(80,'\n');
      fin >> r2;      fin.ignore(80,'\n');
      fin >> r3;      fin.ignore(80,'\n');
      fin >> r4;      fin.ignore(80,'\n');
      fin >> r5;      fin.ignore(80,'\n');
      fin >> r6;      fin.ignore(80,'\n');
                      fin.ignore(80,'\n');
      fin >> uF ;     fin.ignore(80,'\n');
                      fin.ignore(80,'\n');
      fin >> u1;      fin.ignore(80,'\n');
      fin >> u2;      fin.ignore(80,'\n');
      fin >> u3;      fin.ignore(80,'\n');
      fin >> u4;      fin.ignore(80,'\n');
      fin >> u5;      fin.ignore(80,'\n');
      fin >> u6;      fin.ignore(80,'\n');
                      fin.ignore(80,'\n');
      fin >> wF;      fin.ignore(80,'\n');
                      fin.ignore(80,'\n');
      fin >> w1;      fin.ignore(80,'\n');
      fin >> w2;      fin.ignore(80,'\n');
      fin >> w3;      fin.ignore(80,'\n');
      fin >> w4;      fin.ignore(80,'\n');
                      fin.ignore(80,'\n');
      fin >> oF;      fin.ignore(80,'\n');
                      fin.ignore(80,'\n');
      fin >> o1;      fin.ignore(80,'\n');
      fin >> o2;      fin.ignore(80,'\n');
      fin >> o3;      fin.ignore(80,'\n');
      fin >> o4;      fin.ignore(80,'\n');
                      fin.ignore(80,'\n');
      fin >> imax;    fin.ignore(80,'\n');
      fin >> nmax;    fin.ignore(80,'\n');
                      fin.ignore(80,'\n');
      fin >> IMAX;    fin.ignore(80,'\n');
      fin >> delta;   fin.ignore(80,'\n');
                      fin.ignore(80,'\n');
      fin >> mesh;    fin.ignore(80,'\n');
                      fin.ignore(80,'\n');
      fin >> RBW;     fin.ignore(80,'\n');
      fin >> RBO;     fin.ignore(80,'\n');
                      fin.ignore(80,'\n');
      fin >> ED;      fin.ignore(80,'\n');
      fin >> AZ  ;    fin.ignore(80,'\n');

      fin.close();

      if (!(CppVersion==InpVersion))  // Versionskontrolle
     {
      cout << "\n >> Warnung: Unterschiedliche Versionsnummern in \"konti.cpp\" und \"input.dat\" !\n\n";
      abort();
     }

     }

   /* -----------------------------------------------------------------------
      Ende von ParameterEinlesen
      ----------------------------------------------------------------------- */



   /* -----------------------------------------------------------------------
      Anfang von ParameterAusgeben
      ----------------------------------------------------------------------- */

      void ParameterAusgeben(void)
     {

      cout << setiosflags(ios::left)

           << " xa   = " << setw(15) << xa
           << " xe   = " << setw(15) << xe
           << " dx   = " << setw(15) << dx[0]
           << " imax = "             << imax << "\n\n"

           << " ta   = " << setw(15) << ta
           << " te   = " << setw(15) << te
           << " dt   = " << setw(15) << dt
           << " nmax = " << setw(15) << nmax << "\n\n"

           << " r1   = " << setw(15) << r1
           << " r2   = " << setw(15) << r2
           << " r3   = " << setw(15) << r3
           << " r4   = " << setw(15) << r4
           << " r5   = " << setw(15) << r5
           << " r6   = "             << r6   << "\n\n"

           << " u1   = " << setw(15) << u1
           << " u2   = " << setw(15) << u2
           << " u3   = " << setw(15) << u3
           << " u4   = " << setw(15) << u4
           << " u5   = " << setw(15) << u5
           << " u6   = "             << u6   << "\n\n"

           << " w1   = " << setw(15) << w1
           << " w2   = " << setw(15) << w2
           << " w3   = " << setw(15) << w3
           << " w4   = "             << w4   << "\n\n"

           << " o1   = " << setw(15) << o1
           << " o2   = " << setw(15) << o2
           << " o3   = " << setw(15) << o3
           << " o4   = "             << o4   << "\n\n"

           << " dx/dt= " << setw(15) << dx[0]/dt

           << " t0   = " << setw(15) << t0
           << " ED   = " << setw(15) << ED
           << " mesh = " << setw(15) << mesh << "\n\n"

           << " IMAX = " << setw(15) << IMAX
           << " delta= " << setw(15) << delta << "\n\n"

           << " AZ   = " << setw(15) << AZ   << "\n\n";

   /* Randbedingungen und Funktionsnamen ausgeben */

      RandbedingungAusgeben(RBW,"RBW");

      RandbedingungAusgeben(RBO,"RBO");
                                         if ( RBW == DynamicBoundary )
      FunktionsNameAusgeben(wF ,"wF" );
                                         if ( RBO == DynamicBoundary )
      FunktionsNameAusgeben(oF ,"oF" );
                                     
      FunktionsNameAusgeben(rF ,"rF" );
                                     
      FunktionsNameAusgeben(uF ,"uF" );

   /* Warnmeldungen */

      if ( RBW == DirichletBoundary && u[imin-1] < 0.0 )
     {
      cout << "  >> Warnung: u(xa) < 0 ! (Reversed Flow)\n\n";
     }

      if ( RBO == DirichletBoundary && u[imax+1] > 0.0 )
     {
      cout << "  >> Warnung: u(xe) > 0 ! (Reversed Flow)\n\n";
     }

      if ( RBW == NeumannBoundary && u[imin-1] < 0.0 )
     {
      cout << "  >> Warnung: u(xa) < 0 ! (Reversed Flow)\n\n";
     }

      if ( RBO == NeumannBoundary && u[imax+1] > 0.0 )
     {
      cout << "  >> Warnung: u(xe) > 0 ! (Reversed Flow)\n\n";
     }

      if ( RBW == OutletBoundary && u[imin] > 0.0 )
     {
      cout << "  >> Warnung: u(xa) > 0 ! (Reversed Flow)\n\n";
     }

      if ( RBO == OutletBoundary && u[imax] < 0.0 )
     {
      cout << "  >> Warnung: u(xe) < 0 ! (Reversed Flow)\n\n";
     }

     }

   /* -----------------------------------------------------------------------
      Ende von ParameterAusgeben
      ----------------------------------------------------------------------- */



   /* -----------------------------------------------------------------------
      Anfang von RandbedingungAusgeben
      ----------------------------------------------------------------------- */

      void RandbedingungAusgeben(int Nummer,char *Text)
     {
      
      switch(Nummer)
     {
       case  0: cout << " " << Text << "   = 0: Wall-Boundary               \n\n";  break;
       case  1: cout << " " << Text << "   = 1: Dirichlet-Boundary          \n\n";  break;
       case  2: cout << " " << Text << "   = 2: Neumann-Boundary            \n\n";  break;
       case  3: cout << " " << Text << "   = 3: Periodic-Boundary           \n\n";  break;
       case  4: cout << " " << Text << "   = 4: Dynamic-Boundary            \n\n";  break;
       case  5: cout << " " << Text << "   = 5: Outlet-Boundary             \n\n";  break;
     }

     }

   /* -----------------------------------------------------------------------
      Ende von RandbedingungAusgeben
      ----------------------------------------------------------------------- */



   /* -----------------------------------------------------------------------
      Anfang von FunktionsNameAusgeben
      ----------------------------------------------------------------------- */

      void FunktionsNameAusgeben(int Nummer,char *Text)
     {
      
      switch(Nummer)
     {
       case  0: cout << " " << Text << "   =  0: Konstantefunktion           \n\n";  break;
       case  1: cout << " " << Text << "   =  1: Lineare Funktion            \n\n";  break;
       case  2: cout << " " << Text << "   =  2: Parabelfunktion             \n\n";  break;
       case  3: cout << " " << Text << "   =  3: Rechteckfunktion            \n\n";  break;
       case  4: cout << " " << Text << "   =  4: Dreieckfunktion             \n\n";  break;
       case  5: cout << " " << Text << "   =  5: Saegezahnfunktion           \n\n";  break;
       case  6: cout << " " << Text << "   =  6: Lineare Rampe               \n\n";  break;
       case  7: cout << " " << Text << "   =  7: Cosinus Rampe               \n\n";  break;
       case  8: cout << " " << Text << "   =  8: Sinusfunktion               \n\n";  break;
       case  9: cout << " " << Text << "   =  9: Cosinusfunktion             \n\n";  break;
       case 10: cout << " " << Text << "   = 10: Exponentialfunktion         \n\n";  break;
       case 11: cout << " " << Text << "   = 11: Gaussfunktion               \n\n";  break;
       case 12: cout << " " << Text << "   = 12: Diracfunktion               \n\n";  break;
       case 13: cout << " " << Text << "   = 12: Heavisidefunktion           \n\n";  break;
     }

     }

   /* -----------------------------------------------------------------------
      Ende von FunktionsNameAusgeben
      ----------------------------------------------------------------------- */



   /* -----------------------------------------------------------------------
      Anfang von SimulationVorbereiten
      ----------------------------------------------------------------------- */

      void SimulationVorbereiten(void)
     {

      int i;          // Lokaler Schleifenzähler

      imin = 1;       // Erster Zellmittelpunkt

      NMAX = 0;       // Zähler auf Null setzen

      ifstream fin;   // Objekt für Dateieingabe
      ofstream fout;  // Objekt für Dateiausgabe

      fout << setiosflags(ios::scientific) << setprecision(13);

   /* ----------------------------------------------------------------------- */

   /* Speicher reservieren */

      x   = new double[imax+2];
     dx   = new double[imax+2];

      u   = new double[imax+2];
      rho = new double[imax+2];

      f   = new double[imax+2];

      aW  = new double[imax];
      aP  = new double[imax];
      aE  = new double[imax];

   /* ----------------------------------------------------------------------- */

   /* Inkonsistente Daten aus input.dat überschreiben */

      if ( RBW == PeriodicBoundary ) RBO = PeriodicBoundary;
      else
      if ( RBO == PeriodicBoundary ) RBW = PeriodicBoundary;

   /* ----------------------------------------------------------------------- */

   /* Diskretisierungsschema */

      if (RBW == PeriodicBoundary  ) { ImpliziterZeitschritt = ImpliziterZeitschrittPBC;
                                       HilfsfelderBerechnen  = HilfsfelderBerechnenPBC ; }

      if (RBW != PeriodicBoundary  ) { ImpliziterZeitschritt = ImpliziterZeitschrittSTD;
                                       HilfsfelderBerechnen  = HilfsfelderBerechnenSTD ; }

   /* ----------------------------------------------------------------------- */

   /* Äquidistantes Gitter berechnen */

      if(mesh==0)
     {
      for (i=imin-1;i<=imax+1;i++) dx[i] = (xe-xa)/imax;    // Gitterabstände
                                    x[0] = xa - dx[0]/2.0;  // Geisterknoten West
      for (i=imin  ;i<=imax+1;i++)  x[i] = x[0] + i*dx[0];  // Zellmittelpunkte
     }

   /* Externes Gitter einlesen */

      if(mesh==1)
     {
      fin.open("mesh.dat");     // Zellmittelpunkte und Gitterabstände
      for (i=imin-1;i<=imax+1;i++) { fin >> x[i]; fin >> dx[i]; }
      fin.close();
     }

   /* ----------------------------------------------------------------------- */

      N  = 1 + nmax/1000;

      Pi = 2.0*acos(0.0);

   /* Diverse Zeitpunkte berechnen */

      if (AZ==0)          // Mit Anfangszustand starten
     {
      t0 = ta;
     }
      else                // Simulation fortsetzen
     {
      fin.open("te.out");
      fin >> t0;          // Globale Anfangszeit t0 einlesen
      fin >> ta;          // Lokale  Anfangszeit ta einlesen
      fin.close();
     }

      te = ta + nmax*dt;  // Lokalen Endzeitpunkt te berechnen

   /* Diverse Zeitpunkte speichern */

      fout.open("te.out");
      fout << t0 << "\n";
      fout << te << "\n";  // te ist nächstes ta
      fout.close();

   /* ----------------------------------------------------------------------- */

   /* Funktionen für den Anfangszustand und die Ränder definieren */

      U0 = Funktion[uF]; R0 = Funktion[rF];
      W0 = Funktion[wF]; O0 = Funktion[oF];

   /* ------------------------------------------------------------------ */

   /* Anfangszustand definieren */

      if (AZ==0)                    // Neue Simulation starten
     {

      for (i=imin-1;i<=imax+1;i++)  // Alle Zellmittelpunkte
     {
        u[i] = u4*U0( ( x[i]-u1 )/u3 ) + u2;
      rho[i] = r4*R0( ( x[i]-r1 )/r3 ) + r2;      
     }

   /* Anfangszustand für Grafik ohne Ränder speichern */

      fout.open("u.out");
      for (i=imin;i<=imax;i++) fout << x[i] << " " << u[i]   << "\n";
      fout.close();

      fout.open("rho0.out");
      for (i=imin;i<=imax;i++) fout << x[i] << " " << rho[i] << "\n";
      fout.close();

   /* Optional Randzellen anpassen */

        u[imin-1] += u5;
        u[imax+1] += u6;
      rho[imin-1] += r5;
      rho[imax-1] += r6;

      if ( RBW == DynamicBoundary ) rho[imin-1] = w4*W0( ( ta - w1 )/w3 ) + w2;
      if ( RBO == DynamicBoundary ) rho[imax+1] = o4*O0( ( te - r1 )/o3 ) + o2;

      if ( RBW == PeriodicBoundary )
     {

        u[imin-1] =   u[imax];
        u[imax+1] =   u[imin];

      rho[imin-1] = rho[imax];
      rho[imax+1] = rho[imin];

     }

     }

   /* ------------------------------------------------------------------ */

      if (AZ==1)  // Simulation fortsetzen: u und rho einlesen
     {

      fin.open("u.out");
      for(i=imin;i<=imax;i++) { fin >> u[i]; fin >> u[i];}      // Ortskoordinate überlesen
      fin.close();

      fin.open("rho.out");
      for(i=imin;i<=imax;i++) { fin >> rho[i]; fin >> rho[i];}   // Ortskoordinate überlesen
      fin.close();

   /* Geisterzellen separat einlesen */

      fin.open("Boundary.out");
      fin >>   u[imin-1];
      fin >>   u[imax+1];
      fin >> rho[imin-1];
      fin >> rho[imax+1];
      fin.close();

     }

   /* ----------------------------------------------------------------------- */

   /* Diverse Größen berechnen*/

      MasseBerechnen();
      EkinBerechnen();
      ImpulsBerechnen();
      HilfsfelderBerechnen();

   /* Bei NeumannBoundary die erste Ableitung speichern */

      if ( RBW == NeumannBoundary ) GetGradient( rho, imin, dxrhoW);  // RBW merken
      if ( RBO == NeumannBoundary ) GetGradient( rho, imax, dxrhoO);  // RBO merken

     }

   /* -----------------------------------------------------------------------
      Ende von SimulationVorbereiten
      ----------------------------------------------------------------------- */



   /* -----------------------------------------------------------------------
      Anfang von SimulationDurchfuehren
      ----------------------------------------------------------------------- */

      void SimulationDurchfuehren(void)
     {

   /* ------------------------------------------------------------------ */

   /* Lokale Variablen deklarieren */

      int i ;                      // Schleifenzähler
                                   
      unsigned long n;             // Schleifenzähler

      ofstream Mout, Ekin, pxOut;  // Objekte für die Dateiausgabe

   /* ------------------------------------------------------------------ */

      cout << setiosflags(ios::left);

      Mout  << setiosflags(ios::scientific) << setprecision(13);
      Ekin  << setiosflags(ios::scientific) << setprecision(13);
      pxOut << setiosflags(ios::scientific) << setprecision(13);
      Dout << setiosflags(ios::scientific) << setprecision(13);

   /* ------------------------------------------------------------------ */

   /* Simulation durchführen */

      if (AZ==1)
     {

       Mout.open("M.out",   ios::app);
       Ekin.open("Ekin.out",ios::app);
      pxOut.open("px.out",  ios::app);
       Dout.open("D.out",ios::app);

     }
      else
     {

       Mout.open("M.out"   );
       Ekin.open("Ekin.out");
      pxOut.open("px.out"  );
       Dout.open("D.out");

       Mout << ta << " " << M  << "\n";
       Ekin << ta << " " << Ek << "\n";
      pxOut << ta << " " << px << "\n";

     }

      cout << " Status: 0%\r" << flush;

   /* Zeitschleife */

      for (n=1;n<=nmax;n++)
     {

      t = ta + n*dt;

   /* Euler-Zeitschritt */

      ImpliziterZeitschritt();

   /* ----------------------------------------------------------------------- */

      if(n%N==0)
     {
      MasseBerechnen();
      EkinBerechnen();
      ImpulsBerechnen();
       Mout << t << " " << M  << "\n";
       Ekin << t << " " << Ek << "\n";
      pxOut << t << " " << px << "\n";
      cout << " Status: " << n*100/nmax << "%\r" << flush;
     }

     }  // Ende Zeitschleife

   /* ----------------------------------------------------------------------- */

      if( nmax%N != 0 || nmax == 0)
     {
      MasseBerechnen();
      EkinBerechnen();
      ImpulsBerechnen();
       Mout << te << " " << M  << "\n";
       Ekin << te << " " << Ek                << "\n";
      pxOut << te << " " << px << "\n";
     }

       Mout.close();
       Ekin.close();
      pxOut.close();
       Dout.close();

      cout << " Status: 100%\n\n" << flush;

      cout << " IMAX " << NMAX << " mal erreicht!\n\n";

     }

   /* -----------------------------------------------------------------------
      Ende von SimulationDurchfuehren
      ----------------------------------------------------------------------- */



   /* -----------------------------------------------------------------------
      Anfang von ErgebnisseSpeichern
      ----------------------------------------------------------------------- */

      void ErgebnisseSpeichern(void)
     {

      int i;  // Schleifenzähler

      ofstream fout;

      fout << setiosflags(ios::scientific) << setprecision(13);

      fout.open("C.out");       // Convective Courant Number
      for (i=imin;i<=imax;i++) fout << x[i] << " " << fabs(u[i])*dt/dx[i] << "\n";
      fout.close();

      fout.open("rho.out");
      for (i=imin;i<=imax;i++) fout << x[i] << " " << rho[i] << "\n";
      fout.close();

      fout.open("jm.out");
      for (i=imin;i<=imax;i++) fout << x[i] << " " << rho[i]*u[i] << "\n";
      fout.close();

   /* Geisterzellen separat speichern */

      fout.open("Boundary.out");
      fout <<   u[imin-1] << "\n";
      fout <<   u[imax+1] << "\n";
      fout << rho[imin-1] << "\n";
      fout << rho[imax+1] << "\n";
      fout.close();

     }

   /* -----------------------------------------------------------------------
      Ende von ErgebnisseSpeichern
      ----------------------------------------------------------------------- */



   /* -----------------------------------------------------------------------
      Anfang von ErgebnisseDarstellen
      ----------------------------------------------------------------------- */

      void ErgebnisseDarstellen(void)
     {

      if (ED==1) system("showresults");

     }

   /* -----------------------------------------------------------------------
      Ende von ErgebnisseDarstellen
      ----------------------------------------------------------------------- */

   
   
   /* -----------------------------------------------------------------------
      Anfang von ProgrammBeenden
      ----------------------------------------------------------------------- */

      void ProgrammBeenden(void)
     {

      cout << " Programm beendet.\n\n";

     }

   /* -----------------------------------------------------------------------
      Ende von ProgrammBeenden
      ----------------------------------------------------------------------- */



   /* ------------------------------------------------------------------
      Anfang von VolumenIntegration
      ------------------------------------------------------------------ */

      double VolumenIntegration(double *m)  // Rechteckregel (zellzentriert)
     {

      int i;      // Lokaler Schleifenzähler

      double M;   // Lokale Extensive Größe M

      M = 0.0;

      for (i=imin;i<=imax;i++) M += m[i]*dx[i];

      return M;

     }

   /* ------------------------------------------------------------------
      Ende von VolumenIntegration
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von HilfsfelderBerechnenPBC
      ------------------------------------------------------------------ */

      void HilfsfelderBerechnenPBC(void)
     {

      int i;                    // Lokaler Schleifenzähler

      for (i=imin;i<=imax;i++)  // Standard für alle Zellen
     {
      aW[i] = - ( u[i-1] + fabs(u[i-1]) ) * dt/dx[i]/2.0 ;
      aP[i] =       1.0  + fabs(u[ i ])   * dt/dx[i]     ;
      aE[i] =   ( u[i+1] - fabs(u[i+1]) ) * dt/dx[i]/2.0 ;
     }


     }

   /* ------------------------------------------------------------------
      Ende von HilfsfelderBerechnenPBC
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von HilfsfelderBerechnenSTD
      ------------------------------------------------------------------ */

      void HilfsfelderBerechnenSTD(void)
     {

      int i;                        // Lokaler Schleifenzähler

      for (i=imin+1;i<=imax-1;i++)  // Standard für innere Zellen
     {
      aW[i] = - ( u[i-1] + fabs(u[i-1]) ) * dt/dx[i]/2.0 ;
      aP[i] =       1.0  + fabs(u[ i ])   * dt/dx[i]     ;
      aE[i] =   ( u[i+1] - fabs(u[i+1]) ) * dt/dx[i]/2.0 ;
     }

   /* ----------------------------------------------------------------------- */

      if ( RBW == WallBoundary )  // Kein Fluss zwischen W und P
     {

         i  = imin                                            ;
      aW[i] = 0.0                                             ;
      aP[i] = 1.0  + ( u[ i ] + fabs(u[ i ]) ) * dt/dx[i]/2.0 ;
      aE[i] =        ( u[i+1] - fabs(u[i+1]) ) * dt/dx[i]/2.0 ;
      
     }

      if ( RBO == WallBoundary )  // Kein Fluss zwischen E und P
     {

         i  = imax                                            ;
      aW[i] =      - ( u[i-1] + fabs(u[i-1]) ) * dt/dx[i]/2.0 ;
      aP[i] = 1.0  - ( u[ i ] - fabs(u[ i ]) ) * dt/dx[i]/2.0 ;
      aE[i] = 0.0  ;
      
     }

   /* ----------------------------------------------------------------------- */

      if ( RBW == DirichletBoundary || RBW == DynamicBoundary || RBW == NeumannBoundary )  // Kein Fluss von P nach W
     {

         i  = imin                                            ;
      aW[i] =      - ( u[i-1] + fabs(u[i-1]) ) * dt/dx[i]/2.0 ;
      aP[i] = 1.0  + ( u[ i ] + fabs(u[ i ]) ) * dt/dx[i]/2.0 ;
      aE[i] =        ( u[i+1] - fabs(u[i+1]) ) * dt/dx[i]/2.0 ;

     }

      if ( RBO == DirichletBoundary || RBO == DynamicBoundary || RBO == NeumannBoundary )  // Kein Fluss von P nach E
     {

         i  = imax                                            ;
      aW[i] =      - ( u[i-1] + fabs(u[i-1]) ) * dt/dx[i]/2.0 ;
      aP[i] = 1.0  - ( u[ i ] - fabs(u[ i ]) ) * dt/dx[i]/2.0 ;
      aE[i] =        ( u[i+1] - fabs(u[i+1]) ) * dt/dx[i]/2.0 ;
      
     }

   /* ----------------------------------------------------------------------- */

      if ( RBW == OutletBoundary )     // Kein Fluss von W nach P
     {

          i = imin                                       ;
      aW[i] =       0.0                                  ;
      aP[i] =       1.0  + fabs(u[ i ])   * dt/dx[i]     ;
      aE[i] =   ( u[i+1] - fabs(u[i+1]) ) * dt/dx[i]/2.0 ;

     }

      if ( RBO == OutletBoundary )     // Kein Fluss von E nach P
     {
          i = imax                                       ;
      aW[i] = - ( u[i-1] + fabs(u[i-1]) ) * dt/dx[i]/2.0 ;
      aP[i] =       1.0  + fabs(u[ i ])   * dt/dx[i]     ;
      aE[i] =       0.0                                  ;

     }

   /* ----------------------------------------------------------------------- */

     }

   /* ------------------------------------------------------------------
      Ende von HilfsfelderBerechnenSTD
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von MasseBerechnen
      ------------------------------------------------------------------ */

      void MasseBerechnen(void)
     {

      M = VolumenIntegration(rho);

     }

   /* ------------------------------------------------------------------
      Ende von MasseBerechnen
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von EkinBerechnen
      ------------------------------------------------------------------ */

      void EkinBerechnen(void)     // Bezogen auf CellCenter!
     {

      int i;  // Schleifenzähler

      for (i=imin;i<=imax;i++) f[i] = 0.5*rho[i]*u[i]*u[i];

      Ek = VolumenIntegration(f);

     }

   /* ------------------------------------------------------------------
      Ende von EkinBerechnen
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von ImpulsBerechnen
      ------------------------------------------------------------------ */

      void ImpulsBerechnen(void)     // Bezogen auf CellCenter!
     {

      int i;  // Schleifenzähler

      for (i=imin;i<=imax;i++) f[i] = rho[i]*u[i];

      px = VolumenIntegration(f);

     }

   /* ------------------------------------------------------------------
      Ende von ImpulsBerechnen
      ------------------------------------------------------------------ */



   /* -----------------------------------------------------------------------
      Anfang von ImpliziterZeitschrittPBC
      ----------------------------------------------------------------------- */

      void ImpliziterZeitschrittPBC(void)
     {

      int i, K ; // Schleifenzähler

      double D, df;

   /* Startlösung zu K = 0 für Gauß-Seidel erzeugen */

      for (i=imin-1;i<=imax+1;i++)
     {
      f[i] = rho[i];
     }

   /* ------------------------------------------------------------------ */

      for (K=0;K<IMAX;K++)          // Anfang innere Iterationen
     {                              
                                    
      D  = 0.0;                     // Defekt auf Null setzen
                                    
      f[imin-1] = f[imax];          // Geisterzelle West ohne Defektberechnung (Pragmatismus)
                                    
      for (i=imin;i<=imax;i++)      // Innere Zellen mit Defektberechnung
     {
      df = -( aW[i]*f[i-1] + aE[i]*f[i+1] - rho[i]   ) /aP[i] - f[i];
      f[i] += df;
      D  = (fabs(df)>D) * fabs(df) + (fabs(df)<=D) * D;   // if-Verzweigung vermeiden
     }

      f[imax+1] = f[imin];          // Geisterzelle Ost ohne Defektberechnung (Pragmatismus)
                                    
      if ( D < delta ) break;       // D = 0 in Ausgabe vermeiden, daher hier break!
                                    
      Dout << D << "\n";            // ln-Skala: D > 0 hier!
                                    
     }                              // Ende innere Iterationen

   /* ------------------------------------------------------------------ */

      for (i=imin-1;i<=imax+1;i++)  // Rückspeichern
     {
      rho[i] = f[i];
     }

      if (K==IMAX) NMAX++ ;

     }

   /* -----------------------------------------------------------------------
      Ende von ImpliziterZeitschrittPBC
      ----------------------------------------------------------------------- */



   /* -----------------------------------------------------------------------
      Anfang von ImpliziterZeitschrittSTD
      ----------------------------------------------------------------------- */

      void ImpliziterZeitschrittSTD(void)
     {

      int i, K ; // Schleifenzähler

      double D, df;

   /* Startlösung für Gauß-Seidel erzeugen */

      for (i=imin-1;i<=imax+1;i++)
     {
      f[i] = rho[i];
     }

   /* DynamicBoundary für neuen Zeitschritt auswerten */

      if ( RBW == DynamicBoundary ) f[imin-1] = w4*W0( ( t - w1 )/w3 ) + w2;
      if ( RBO == DynamicBoundary ) f[imax+1] = o4*O0( ( t - o1 )/o3 ) + o2;

      for (K=0;K<IMAX;K++)  // Innere Iterationen
     {

      D  = 0.0;  // Geisterzellen ohne Defektberechnung (Pragmatismus)

      if ( RBW == NeumannBoundary ) SetGradient( f, imin, dxrhoW );

   /* Die i-te Gleichung: aW[i]*f[i-1] + aP[i]*f[i] + aE[i]*f[i+1] = rho[i] */

      for (i=imin;i<=imax;i++)
     {
      df = ( rho[i] - aW[i]*f[i-1] - aE[i]*f[i+1] ) / aP[i] - f[i];
      f[i] += df;
      if (fabs(df)>D ) D  = fabs(df);
     }

      if ( RBO == NeumannBoundary ) SetGradient( rho, imax, dxrhoO );

      if ( D < delta ) break;

      Dout << D << "\n";  // ln-Skala: D = 0 vermeiden, daher Dout nach break

     }

      for (i=0;i<=imax;i++)  // Rückspeichern
     {
      rho[i] = f[i];
     }

      if (K==IMAX) NMAX++ ;

     }

   /* -----------------------------------------------------------------------
      Ende von ImpliziterZeitschrittSTD
      ----------------------------------------------------------------------- */



   /* ------------------------------------------------------------------
      Anfang von Konstante
      ------------------------------------------------------------------ */

      double Konstante(double x)
     {
      return 1.0;
     }

   /* ------------------------------------------------------------------
      Ende von Konstante
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von Linear
      ------------------------------------------------------------------ */

      double Linear(double x)
     {
      return x;
     }

   /* ------------------------------------------------------------------
      Ende von Linear
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von Parabel
      ------------------------------------------------------------------ */

      double Parabel(double x)
     {
      return x*x;
     }

   /* ------------------------------------------------------------------
      Ende von Parabel
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von Rechteck
      ------------------------------------------------------------------ */

      double Rechteck(double x)
     {
      if ( x <= -0.5 ) return 0.0;  // Wegen Normierung auf 1 hier <=
      if ( x >   0.5 ) return 0.0;
                       return 1.0;
     }

   /* ------------------------------------------------------------------
      Ende von Rechteck
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von Dreieck
      ------------------------------------------------------------------ */

      double Dreieck(double x)
     {
      if ( x <= -1.0 ) return 0.0;
      if ( x <=  0.0 ) return 1.0 + x;
      if ( x <= +1.0 ) return 1.0 - x;
                       return 0.0;
     }

   /* ------------------------------------------------------------------
      Ende von Dreieck
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von Saegezahn
      ------------------------------------------------------------------ */

      double Saegezahn(double x)
     {
      if ( x <=  0.0 ) return 0.0;
      if ( x <= +1.0 ) return x;
                       return 0.0;
     }

   /* ------------------------------------------------------------------
      Ende von Saegezahn
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von Linearrampe
      ------------------------------------------------------------------ */

      double Linearrampe(double x)
     {
      if ( x <=  0.0 ) return 0.0;
      if ( x <= +1.0 ) return x;
                       return 1.0;
     }

   /* ------------------------------------------------------------------
      Ende von Linearrampe
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von Cosinusrampe
      ------------------------------------------------------------------ */

      double Cosinusrampe(double x)
     {
      if ( x <=  0.0 ) return 0.0;
      if ( x <= +1.0 ) return 0.5*(1.0-cos(Pi*x));
                       return 1.0;
     }

   /* ------------------------------------------------------------------
      Ende von Cosinusrampe
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von Sinus
      ------------------------------------------------------------------ */

      double Sinus(double x)
     {
      return sin(2.0*Pi*x);
     }

   /* ------------------------------------------------------------------
      Ende von Sinus
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von Cosinus
      ------------------------------------------------------------------ */

      double Cosinus(double x)
     {
      return cos(2.0*Pi*x);
     }

   /* ------------------------------------------------------------------
      Ende von Cosinus
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von CosPeak
      ------------------------------------------------------------------ */

      double CosPeak(double x)
     {
      if ( x <= -0.5 ) return 0.0;
      if ( x >=  0.5 ) return 0.0;
                       return 0.5*(1.0+cos(2.0*Pi*x));
     }

   /* ------------------------------------------------------------------
      Ende von CosPeak
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von Exponential
      ------------------------------------------------------------------ */

      double Exponential(double x)
     {
      return exp(x);
     }

   /* ------------------------------------------------------------------
      Ende von Exponential
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von Gauss
      ------------------------------------------------------------------ */

      double Gauss(double x) // Funktion ist normiert!
     {

      return exp(-x*x*Pi);

     }

   /* ------------------------------------------------------------------
      Ende von Gauss
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von Dirac
      ------------------------------------------------------------------ */

      double Dirac(double x)           // Funktion ist für äquidistantes
     {                                 //  Gitter normiert!

      if ( x == 0.0 ) return 1.0/dx[0];
                      return 0.0;

     }

   /* ------------------------------------------------------------------
      Ende von Dirac
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von Heaviside
      ------------------------------------------------------------------ */

      double Heaviside(double x)
     {
      if ( x < 0.0 ) return 0.0; 
                     return 1.0;
     }

   /* ------------------------------------------------------------------
      Ende von Heaviside
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von UserDefined01
      ------------------------------------------------------------------ */

      double UserDefined01(double x)
     {
      double y; 
     #include "./UserDef01.txt"   
      ;
      return(y);
     }

   /* ------------------------------------------------------------------
      Ende von UserDefined01
      ------------------------------------------------------------------ */



   /* ------------------------------------------------------------------
      Anfang von UserDefined02
      ------------------------------------------------------------------ */

      double UserDefined02(double x)
     {
      double y; 
     #include "./UserDef02.txt"   
      ;
      return(y);
     }

   /* ------------------------------------------------------------------
      Ende von UserDefined02
      ------------------------------------------------------------------ */



   /* -----------------------------------------------------------------------
      Anfang von GetGradient
      ----------------------------------------------------------------------- */

      void GetGradient(double *m, int i, double &Value)
     {

      if (i == imin) Value =  ( m[ i ] - m[i-1] ) / ( x[ i ] - x[i-1] );
      if (i == imax) Value =  ( m[i+1] - m[ i ] ) / ( x[i+1] - x[ i ] );

     }

   /* -----------------------------------------------------------------------
      Ende von GetGradient
      ----------------------------------------------------------------------- */



   /* -----------------------------------------------------------------------
      Anfang von SetGradient
      ----------------------------------------------------------------------- */

      void SetGradient(double *m, int i, double Value)
     {

      if (i==imin) m[i-1] = m[i] - Value * ( x[ i ] - x[i-1] );
      if (i==imax) m[i+1] = m[i] + Value * ( x[i+1] - x[ i ] );

     }

   /* -----------------------------------------------------------------------
      Ende von SetGradient
      ----------------------------------------------------------------------- */
