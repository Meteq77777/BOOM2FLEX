#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
WRF-FLEXPART-AIRA Tool (Updated Version - Latest Formula v4)
================================================================================
A tool for generating FLEXPART input files for atmospheric dispersion modeling
of explosive events, integrating with WRF meteorological data.

Updated to support multiple TNT calculation methods via WOPTIONS_TYPES.
FIXED: CSK field type changed to FLOAT, range corrected to (0.8~2.0)
FIXED: Crater calculation formula updated: CRATER = 0.4 × TNT^0.408 (CSK no longer used)
FIXED: TNT result converted to kt by dividing by 1,000,000
CORRECTED: Method 3 updated to correct MPV formula: MPV = 150/(R/TNT^0.3)^1.5

Author: Generated for atmospheric modeling applications
License: GPL-3.0
Version: 1.5 (Latest Formula v4 - MPV Formula Corrected)
================================================================================
"""

import os
import sys
import math
import argparse
import configparser
from datetime import datetime
import shutil

class FlexpartAiraProcessor:
    """Main class for processing FLEXPART-AIRA input generation."""
    
    def __init__(self, config_file="namelist.flexpart"):
        """Initialize the processor with configuration file."""
        self.config_file = config_file
        self.config = configparser.ConfigParser()
        self.template_file = "flexinput.txt"
        self.output_file = "flexinput_BOOM.txt"
        self.errors = []
        
    def read_config(self):
        """Read configuration from namelist file."""
        if not os.path.exists(self.config_file):
            raise FileNotFoundError(f"Configuration file {self.config_file} not found!")
            
        self.config.read(self.config_file)
        print(f"Reading configuration from {self.config_file}")
        
    def validate_input(self):
        """Validate all input parameters according to specifications."""
        print("Validating input parameters...")
        
        # Get parameters from config
        params = {}
        for section in self.config.sections():
            for key, value in self.config.items(section):
                params[key] = value
        
        # Validation rules
        try:
            # Time interval validation
            int_val = float(params.get('int', 0))
            if int_val < 3600:
                self.errors.append("INTERVAL ERROR：TOO SHORT")
                
            # Longitude validations
            llx = float(params.get('llx', 0))
            if not (-180 <= llx <= 180):
                self.errors.append("OUT OF RANGE")
                
            outx = float(params.get('outx', 0))
            if not (-180 <= outx <= 180) or outx <= llx:
                self.errors.append("OUT OF RANGE")
                
            # Latitude validations
            lly = float(params.get('lly', 0))
            if not (-90 <= lly <= 90):
                self.errors.append("OUT OF RANGE")
                
            outy = float(params.get('outy', 0))
            if not (-90 <= outy <= 90) or outy <= lly:
                self.errors.append("OUT OF RANGE")
                
            # Grid point validations
            numx = int(params.get('numx', 0))
            if not (0 <= numx <= 300):
                self.errors.append("OUT OF RANGE")
                
            numy = int(params.get('numy', 0))
            if not (0 <= numy <= 300):
                self.errors.append("OUT OF RANGE")
                
            # Species validation
            specnum = int(params.get('specnum', 1))
            if specnum > 1:
                self.errors.append("SPECIES NUMBER INPUT ERROR, CHANGE 1")
                
            # Check species file exists
            specnam = params.get('specnam', '')
            species_file = f"SPECIES_{specnam}.txt"
            if not os.path.exists(species_file):
                self.errors.append("SPECIES DONT EXIT")
                
            # Release point validation
            numpoint = int(params.get('numpoint', 1))
            if numpoint != 1:
                self.errors.append("RELEASE POINT OUT OF RANGE")
                
            # Central coordinates validation
            cx = float(params.get('cx', 0))
            if not (-180 <= cx <= 180) or not (llx < cx < outx):
                self.errors.append("OUT OF RANGE")
                
            cy = float(params.get('cy', 0))
            if not (-90 <= cy <= 90) or not (lly < cy < outy):
                self.errors.append("OUT OF RANGE")
                
            # WOPTIONS_TYPES validation and dependent parameter validation
            woptions_types = int(params.get('woptions_types', 1))
            if woptions_types not in [1, 2, 3]:
                self.errors.append("TYPES OUT OF RANGE,CHOOSE 1(WTNT) |2(CRATER paraments) |3(MPV)")
            
            # Validate parameters based on WOPTIONS_TYPES
            if woptions_types == 1:
                # Direct TNT input
                tnt = float(params.get('tnt', 0))
                if tnt <= 0:
                    self.errors.append("TNT value must be greater than 0")
                    
            elif woptions_types == 2:
                # Crater parameters - UPDATED: CSK parameter now optional (not used in new formula)
                crater = float(params.get('crater', 0))
                if crater <= 0:
                    self.errors.append("CRATER radius must be greater than 0")
                
                # CSK parameter validation (optional, for backward compatibility)
                csk_raw = params.get('csk', None)
                if csk_raw:
                    # Handle inline comments in config values
                    csk_clean = csk_raw.split('#')[0].strip() if '#' in csk_raw else csk_raw.strip()
                    try:
                        csk = float(csk_clean)
                        if not (0.8 <= csk <= 2.0):
                            self.errors.append("OUT OF RANGE: 1.5-2.0 for HARD SOIL and 0.8-1.2 for ROCK")
                        else:
                            # Provide guidance based on csk value
                            if 1.5 <= csk <= 2.0:
                                print(f"CSK={csk} - suitable for HARD SOIL (Note: CSK not used in current formula)")
                            elif 0.8 <= csk <= 1.2:
                                print(f"CSK={csk} - suitable for ROCK (Note: CSK not used in current formula)")
                            else:
                                print(f"CSK={csk} - intermediate value (Note: CSK not used in current formula)")
                    except ValueError:
                        self.errors.append("CSK must be a valid float number")
                
            elif woptions_types == 3:
                # Maximum Particle Velocity (MPV) - CORRECTED according to image
                mpv = float(params.get('mpv', 0))
                if mpv <= 0:
                    self.errors.append("MPV (Maximum particle velocity) must be greater than 0")
                
                r = float(params.get('r', 0))
                if r <= 0:
                    self.errors.append("R (distance from explosion source) must be greater than 0")
                
            # Particle number validation
            tnr = int(params.get('tnr', 0))
            if not (0 <= tnr <= 50000):
                self.errors.append("TNR must be between 0 and 50000")
                
            # Mass emission validation
            massemit = float(params.get('massemit', 0))
            if massemit <= 0:
                self.errors.append("Mass emission must be greater than 0")
                
            # Time validation
            start_time = params.get('starttimes', '')
            end_time = params.get('endtimes', '')
            id1it1 = params.get('id1it1', '')
            id2it2 = params.get('id2it2', '')
            
            # Convert time strings to datetime objects for comparison
            try:
                start_dt = datetime.strptime(start_time, '%Y%m%d %H%M%S')
                end_dt = datetime.strptime(end_time, '%Y%m%d %H%M%S')
                release_start_dt = datetime.strptime(id1it1, '%Y%m%d %H%M%S')
                release_end_dt = datetime.strptime(id2it2, '%Y%m%d %H%M%S')
                
                if release_start_dt < start_dt:
                    self.errors.append("Release start time must be >= simulation start time")
                if release_end_dt > end_dt:
                    self.errors.append("Release end time must be <= simulation end time")
                    
            except ValueError as e:
                self.errors.append(f"Time format error: {e}")
                
        except (ValueError, KeyError) as e:
            self.errors.append(f"Parameter error: {e}")
            
        # Stop if there are errors
        if self.errors:
            for error in self.errors:
                print(f"ERROR: {error}")
            sys.exit(1)
            
        print("All input parameters validated successfully!")
        return params
    
    def calculate_tnt(self, params):
        """Calculate TNT yield based on WOPTIONS_TYPES."""
        woptions_types = int(params.get('woptions_types', 1))
        
        if woptions_types == 1:
            # Direct TNT input
            tnt = float(params.get('tnt', 0))
            print(f"=== METHOD 1: Direct TNT Input ===")
            print(f"Input TNT value: {tnt} kt")
            print(f"Final calculated TNT: {tnt} kt")
            return tnt
        
        elif woptions_types == 2:
            # Calculate TNT from crater parameters - FIXED: Use latest formula from image
            crater = float(params.get('crater', 0))
            
            print(f"=== METHOD 2: Crater Parameters ===")
            print(f"Input parameters:")
            print(f"  CRATER radius: {crater} m")
            print(f"Calculation steps:")
            print(f"  Formula: CRATER = 0.4 × TNT^0.408")
            print(f"  Solving for TNT: TNT^0.408 = CRATER / 0.4")
            print(f"  TNT = (CRATER / 0.4)^(1/0.408)")
            print(f"  TNT = ({crater} / 0.4)^(1/0.408)")
            
            # FIXED: Use correct formula from latest provided image
            # Formula: CRATER = 0.4 × TNT^0.408
            # Solving for TNT: TNT = (CRATER / 0.4)^(1/0.408)
            intermediate = crater / 0.4
            exponent = 1 / 0.408  # ≈ 2.451
            print(f"  TNT = ({intermediate})^{exponent:.3f}")
            tnt_raw = intermediate ** exponent
            print(f"  TNT (raw calculation) = {tnt_raw:.2e}")
            
            # FIXED: Convert to kt by dividing by 1,000,000
            tnt = tnt_raw / 1000000
            print(f"  TNT (converted to kt) = {tnt_raw:.2e} / 1,000,000 = {tnt:.6f} kt")
            
            print(f"Final calculated TNT: {tnt:.6f} kt")
            print(f"TNT in scientific notation: {tnt:.2e} kt")
            
            # Note: CSK parameter is no longer used in this formula
            csk_raw = params.get('csk', None)
            if csk_raw:
                print(f"Note: CSK parameter ({csk_raw}) present in config but not used in this formula")
            
            return tnt
        
        elif woptions_types == 3:
            # Calculate TNT from Maximum Particle Velocity (MPV) - CORRECTED according to image
            mpv = float(params.get('mpv', 0))  # um/sec
            r = float(params.get('r', 0))      # ft (feet)
            
            print(f"=== METHOD 3: Maximum Particle Velocity (MPV) ===")
            print(f"Input parameters:")
            print(f"  MPV (Maximum particle velocity): {mpv} um/sec")
            print(f"  R (Distance from explosion source): {r} ft")
            print(f"Calculation steps:")
            print(f"  Formula: MPV = 150/(R/TNT^0.3)^1.5")
            print(f"  Solving for TNT:")
            
            # CORRECTED Formula: MPV = 150/(R/TNT^0.3)^1.5
            # Step by step solution:
            # 1. (R/TNT^0.3)^1.5 = 150/MPV
            # 2. R/TNT^0.3 = (150/MPV)^(1/1.5) = (150/MPV)^(2/3)
            # 3. TNT^0.3 = R / (150/MPV)^(2/3)
            # 4. TNT = [R / (150/MPV)^(2/3)]^(1/0.3) = [R / (150/MPV)^(2/3)]^(10/3)
            
            intermediate1 = 150 / mpv
            print(f"  Step 1: 150 / MPV = 150 / {mpv} = {intermediate1:.6f}")
            
            exponent1 = 2/3
            intermediate2 = intermediate1 ** exponent1
            print(f"  Step 2: (150/MPV)^(2/3) = {intermediate1:.6f}^{exponent1:.3f} = {intermediate2:.6f}")
            
            intermediate3 = r / intermediate2
            print(f"  Step 3: R / (150/MPV)^(2/3) = {r} / {intermediate2:.6f} = {intermediate3:.6f}")
            
            exponent2 = 10/3  # 1/0.3 = 10/3
            tnt = intermediate3 ** exponent2
            print(f"  Step 4: TNT = ({intermediate3:.6f})^(10/3) = {intermediate3:.6f}^{exponent2:.3f} = {tnt:.6f}")
            
            print(f"Final calculated TNT: {tnt:.6f} kt")
            print(f"TNT in scientific notation: {tnt:.2e} kt")
            return tnt
        
        else:
            raise ValueError(f"Unsupported WOPTIONS_TYPES: {woptions_types}")
    
    def calculate_heights(self, tnt):
        """Calculate ZL and ZT based on TNT yield."""
        print(f"\n=== CALCULATING RELEASE HEIGHTS ===")
        print(f"Input TNT yield: {tnt:.6f} kt")
        
        # ZL calculation
        if tnt <= 4.07:
            a, b = 2228, 0.3463
            print(f"TNT <= 4.07 kt, using: a={a}, b={b}")
        else:
            a, b = 2661, 0.2198
            print(f"TNT > 4.07 kt, using: a={a}, b={b}")
        zl = a * (tnt ** b)
        print(f"ZL = {a} × {tnt:.6f}^{b} = {zl:.2f} m")
        
        # ZT calculation
        if tnt < 2.29:
            c, d = 3597, 0.2533
            print(f"TNT < 2.29 kt, using: c={c}, d={d}")
        elif 2.29 <= tnt < 19:
            c, d = 3170, 0.4077
            print(f"2.29 <= TNT < 19 kt, using: c={c}, d={d}")
        else:
            c, d = 6474, 0.1650
            print(f"TNT >= 19 kt, using: c={c}, d={d}")
        zt = c * (tnt ** d)
        print(f"ZT = {c} × {tnt:.6f}^{d} = {zt:.2f} m")
        
        print(f"*** Final heights: ZL={zl:.2f} m, ZT={zt:.2f} m ***")
        return zl, zt
    
    def calculate_release_area(self, cx, cy, tnt):
        """Calculate release area coordinates based on TNT yield."""
        print(f"\n=== CALCULATING RELEASE AREA ===")
        print(f"Center coordinates: ({cx}, {cy})")
        print(f"TNT yield: {tnt:.6f} kt")
        
        # Calculate R using the correct formula
        lg_tnt = math.log10(tnt)
        print(f"log10({tnt:.6f}) = {lg_tnt:.6f}")
        
        exponent = 6.7553 + 0.32055 * lg_tnt + 0.01137478 * (lg_tnt ** 2)
        print(f"Exponent = 6.7553 + 0.32055×{lg_tnt:.6f} + 0.01137478×{lg_tnt:.6f}² = {exponent:.6f}")
        
        r_meters = math.exp(exponent)
        print(f"R (meters) = exp({exponent:.6f}) = {r_meters:.2f} m")
        
        # FIXED: Convert R from meters to degrees by dividing by 111000
        r_degrees = r_meters / 111000
        print(f"R (degrees) = {r_meters:.2f} / 111000 = {r_degrees:.6f} degrees")
        
        # Calculate corners using the corrected formula
        # RLLX = CX - R/111000, RLLY = CY - R/111000
        # RTRX = CX + R/111000, RTRY = CY + R/111000
        rllx = cx - r_degrees
        rlly = cy - r_degrees
        rtrx = cx + r_degrees
        rtry = cy + r_degrees
        
        print(f"Release area corners:")
        print(f"  RLLX = CX - R/111000 = {cx} - {r_degrees:.6f} = {rllx:.6f}")
        print(f"  RLLY = CY - R/111000 = {cy} - {r_degrees:.6f} = {rlly:.6f}")
        print(f"  RTRX = CX + R/111000 = {cx} + {r_degrees:.6f} = {rtrx:.6f}")
        print(f"  RTRY = CY + R/111000 = {cy} + {r_degrees:.6f} = {rtry:.6f}")
        print(f"*** Release radius: {r_meters:.2f} meters ({r_degrees:.6f} degrees) ***")
        
        return rllx, rlly, rtrx, rtry
    
    def read_species_data(self, specnam):
        """Read species data from SPECIES_*.txt file."""
        species_file = f"SPECIES_{specnam}.txt"
        print(f"Reading species data from {species_file}")
        
        with open(species_file, 'r') as f:
            lines = f.readlines()
            
        # Find the data line (not header)
        species_data = ""
        for line in lines:
            if not line.startswith('XXXX') and line.strip():
                species_data = line.strip()
                break
                
        return species_data
    
    def generate_flexinput(self, params):
        """Generate the flexinput_BOOM.txt file."""
        print("Generating flexinput_BOOM.txt...")
        
        # Calculate TNT based on method selection
        tnt = self.calculate_tnt(params)
        
        # Get center coordinates
        cx = float(params['cx'])
        cy = float(params['cy'])
        
        # Calculate derived parameters
        zl, zt = self.calculate_heights(tnt)
        rllx, rlly, rtrx, rtry = self.calculate_release_area(cx, cy, tnt)
        
        # Read species data
        species_data = self.read_species_data(params['specnam'])
        
        # Read template file
        if not os.path.exists(self.template_file):
            raise FileNotFoundError(f"Template file {self.template_file} not found!")
            
        with open(self.template_file, 'r') as f:
            template_content = f.read()
        
        # Replace variables in template
        replacements = {
            '$PATH1': params.get('path1', './input/'),
            '$PATH2': params.get('path2', './wrf_data/'),
            '$PATH3': params.get('path3', './available/'),
            '$STARTTIMES': params['starttimes'],
            '$ENDTIMES': params['endtimes'],
            '$INT': params['int'],
            '$LLX': params['llx'],
            '$LLY': params['lly'],
            '$NUMX': params['numx'],
            '$NUMY': params['numy'],
            '$OUTX': params['outx'],
            '$OUTY': params['outy'],
            '$SPECNUM': params['specnum'],
            'SPECNAME': species_data,
            '$NUMPOINT': params['numpoint'],
            '$ID1IT1': params['id1it1'],
            '$ID2IT2': params['id2it2'],
            'RLLX': f"{rllx:.6f}",
            'RLLY': f"{rlly:.6f}",
            'RTRX': f"{rtrx:.6f}",
            'RTRY': f"{rtry:.6f}",
            'ZL': f"{zl:.2f}",
            'ZT': f"{zt:.2f}",
            '$TNR': params['tnr'],
            '$MASSEMIT': params['massemit']
        }
        
        # Apply replacements
        output_content = template_content
        for old, new in replacements.items():
            output_content = output_content.replace(old, str(new))
        
        # Write output file
        with open(self.output_file, 'w') as f:
            f.write(output_content)
        
        print(f"Successfully generated {self.output_file}")
        
        # Print summary
        print("\n" + "="*60)
        print("GENERATION SUMMARY")
        print("="*60)
        print(f"Template file: {self.template_file}")
        print(f"Output file: {self.output_file}")
        print(f"TNT calculation method: {params.get('woptions_types', 1)}")
        print(f"*** FINAL TNT YIELD: {tnt:.6f} kt ***")
        print(f"*** TNT in scientific notation: {tnt:.2e} kt ***")
        print(f"Release heights: ZL={zl:.2f} m, ZT={zt:.2f} m")
        print(f"Release area: ({rllx:.6f}, {rlly:.6f}) to ({rtrx:.6f}, {rtry:.6f})")
        print(f"Release area formula: R/111000 conversion used")
        print(f"Simulation time: {params['starttimes']} to {params['endtimes']}")
        print(f"Species: {params['specnam']}")
        print(f"Particles: {params['tnr']}")
        print(f"Mass emission: {params['massemit']}")
        
        if params.get('woptions_types') == '2':
            print(f"Crater radius: {params.get('crater', 'N/A')} m")
            print(f"Formula used: CRATER = 0.4 × TNT^0.408")
            csk_raw = params.get('csk', None)
            if csk_raw:
                csk_clean = csk_raw.split('#')[0].strip() if '#' in csk_raw else csk_raw.strip()
                print(f"CSK coefficient: {csk_clean} (present but not used in current formula)")
        elif params.get('woptions_types') == '3':
            print(f"MPV (Maximum particle velocity): {params.get('mpv', 'N/A')} um/sec")
            print(f"R (Distance from explosion): {params.get('r', 'N/A')} ft")
            print(f"Formula used: MPV = 150/(R/TNT^0.3)^1.5")
            
        print("="*60)
    
    def run(self):
        """Main execution function."""
        print("="*60)
        print("WRF-FLEXPART-AIRA Tool v1.5 (Latest Formula v4 - MPV Formula Corrected)")
        print("="*60)
        
        try:
            # Read configuration
            self.read_config()
            
            # Validate input parameters
            params = self.validate_input()
            
            # Generate output file
            self.generate_flexinput(params)
            
            print("\nProcessing completed successfully!")
            
        except Exception as e:
            print(f"ERROR: {e}")
            sys.exit(1)

def main():
    """Main function with command line argument parsing."""
    parser = argparse.ArgumentParser(
        description="WRF-FLEXPART-AIRA v1.5 (Latest Formula v4 - MPV Formula Corrected): Generate FLEXPART input files for explosive dispersion modeling",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
TNT Calculation Methods (WOPTIONS_TYPES):
  1 - Direct TNT input (default)
  2 - Calculate from crater parameters (LATEST FORMULA: CRATER = 0.4 × TNT^0.408)
  3 - Calculate from Maximum Particle Velocity (MPV = 150/(R/TNT^0.3)^1.5)

Method 3 Parameters:
  - MPV: Maximum particle velocity (um/sec), must be > 0
  - R: Distance from explosion source (ft), must be > 0

LATEST UPDATES: 
  - Crater calculation formula updated: CRATER = 0.4 × TNT^0.408
  - CSK parameter no longer used in calculation (kept for backward compatibility)
  - TNT result converted to kt by dividing by 1,000,000
  - Release area calculation: R/111000 conversion (R meters to degrees)
  - Method 3 corrected to proper MPV calculation: MPV = 150/(R/TNT^0.3)^1.5
  - Only CRATER radius required for method 2
  - Only MPV and R required for method 3

Release Area Calculation:
  R = exp[6.7553 + 0.32055×Lg(TNT) + 0.01137478×(Lg(TNT))²]
  RLLX = CX - R/111000, RLLY = CY - R/111000
  RTRX = CX + R/111000, RTRY = CY + R/111000

Examples:
  python main.py
  python main.py --config my_namelist.flexpart
  python main.py --template my_template.txt --output my_output.txt
        """
    )
    
    parser.add_argument(
        '--config', '-c',
        default='namelist.flexpart',
        help='Configuration file (default: namelist.flexpart)'
    )
    
    parser.add_argument(
        '--template', '-t',
        default='flexinput.txt',
        help='Template file (default: flexinput.txt)'
    )
    
    parser.add_argument(
        '--output', '-o',
        default='flexinput_BOOM.txt',
        help='Output file (default: flexinput_BOOM.txt)'
    )
    
    parser.add_argument(
        '--version', '-v',
        action='version',
        version='WRF-FLEXPART-AIRA v1.5 (Latest Formula v4 - MPV Formula Corrected)'
    )
    
    args = parser.parse_args()
    
    # Initialize processor
    processor = FlexpartAiraProcessor(args.config)
    processor.template_file = args.template
    processor.output_file = args.output
    
    # Run processing
    processor.run()

if __name__ == "__main__":
    main()