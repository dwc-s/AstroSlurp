"""
Ad Astra - Astronomical Observation Planner

This application allows users to define an observer location and time window,
then queries the SIMBAD astronomical database to find celestial objects
(stars, galaxies, nebulae) that are visible above a specified horizon
and brighter than a specified magnitude limit.
"""
import sys
import csv
import json
import os
import time
import numpy as np
import astropy.units as u
import astroplan
from astroplan import FixedTarget, AltitudeConstraint, AtNightConstraint
from astropy.time import Time
from timezonefinder import TimezoneFinder
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                             QLineEdit, QPushButton, QTextEdit, QFormLayout, QMessageBox,
                             QTableWidget, QTableWidgetItem)

# Mapping of SIMBAD O-types to human-readable descriptions
# Source: https://simbad.cds.unistra.fr/guide/otypes.htx
OTYPE_MAP = {
    '*': 'Star',
    '**': 'Double/Multiple Star',
    'V*': 'Variable Star',
    'PM*': 'High Proper Motion Star',
    'HV*': 'High Velocity Star',
    'Em*': 'Emission-line Star',
    'Be*': 'Be Star',
    'SB*': 'Spectroscopic Binary',
    'EB*': 'Eclipsing Binary',
    'Algol': 'Algol-type Binary',
    'WD*': 'White Dwarf',
    'N*': 'Neutron Star',
    'BH': 'Black Hole',
    'Psr': 'Pulsar',
    'Y*O': 'Young Stellar Object',
    'TT*': 'T Tauri Star',
    'C*': 'Carbon Star',
    'S*': 'S Star',
    'G': 'Galaxy',
    'SyG': 'Seyfert Galaxy',
    'AG': 'Active Galaxy',
    'QSO': 'Quasar',
    'BLLac': 'BL Lac Object',
    'Blazar': 'Blazar',
    'LINER': 'LINER Galaxy',
    'GinCl': 'Galaxy in Cluster',
    'IG': 'Interacting Galaxies',
    'PairG': 'Pair of Galaxies',
    'ClG': 'Cluster of Galaxies',
    'GrG': 'Group of Galaxies',
    'PN': 'Planetary Nebula',
    'Neb': 'Nebula',
    'RNe': 'Reflection Nebula',
    'DkNeb': 'Dark Nebula',
    'MolCld': 'Molecular Cloud',
    'HII': 'HII Region',
    'SNR': 'Supernova Remnant',
    'OpCl': 'Open Cluster',
    'GlCl': 'Globular Cluster',
    'Ass*': 'Association of Stars',
    'Cl*': 'Cluster of Stars',
    'X': 'X-ray Source',
    'IR': 'Infrared Source',
    'Radio': 'Radio Source',
    'UV': 'UV Source',
    'Gam': 'Gamma-ray Source',
    'Pl': 'Planet',
    'Ast': 'Asteroid',
    'Com': 'Comet'
}

class NumericTableWidgetItem(QTableWidgetItem):
    """
    Custom TableItem to ensure numerical sorting for columns like Magnitude and Coordinates.
    Standard QTableWidgetItem sorts alphabetically (e.g., "10" comes before "2").
    """
    def __lt__(self, other):
        try:
            return float(self.text()) < float(other.text())
        except ValueError:
            return super().__lt__(other)

class AdAstraWindow(QMainWindow):
    def __init__(self):
        """
        Initialize the Main Window, Layouts, and Widgets.
        """
        super().__init__()
        self.setWindowTitle("Ad Astra")
        self.resize(1000, 600)

        # Main widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)

        # Left Container
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)

        # Form Layout for inputs
        form_layout = QFormLayout()

        # Config inputs
        self.lat_edit = QLineEdit()
        self.lon_edit = QLineEdit()
        self.elev_edit = QLineEdit()
        self.name_edit = QLineEdit()
        
        form_layout.addRow("Latitude (deg):", self.lat_edit)
        form_layout.addRow("Longitude (deg):", self.lon_edit)
        form_layout.addRow("Elevation (m):", self.elev_edit)
        form_layout.addRow("Place Name:", self.name_edit)

        # User inputs
        self.time_begin_edit = QLineEdit()
        self.time_begin_edit.setPlaceholderText("YYYY-MM-DD HH:MM")
        form_layout.addRow("Start Time:", self.time_begin_edit)

        self.time_end_edit = QLineEdit()
        self.time_end_edit.setPlaceholderText("YYYY-MM-DD HH:MM")
        form_layout.addRow("End Time:", self.time_end_edit)

        self.horizon_edit = QLineEdit()
        form_layout.addRow("Horizon Limit (deg):", self.horizon_edit)

        self.mag_edit = QLineEdit()
        self.mag_edit.setPlaceholderText("e.g. 6.0 (default)")
        form_layout.addRow("Limiting Magnitude:", self.mag_edit)

        self.temp_edit = QLineEdit()
        self.temp_edit.setPlaceholderText("Optional (C)")
        form_layout.addRow("Temperature (C):", self.temp_edit)

        self.press_edit = QLineEdit()
        self.press_edit.setPlaceholderText("Optional (mbar)")
        form_layout.addRow("Pressure (mbar):", self.press_edit)

        self.hum_edit = QLineEdit()
        self.hum_edit.setPlaceholderText("Optional (%)")
        form_layout.addRow("Humidity (%):", self.hum_edit)

        left_layout.addLayout(form_layout)

        # Action Button
        self.create_btn = QPushButton("Create Observer (do this first)")
        self.create_btn.clicked.connect(self.create_observer)
        left_layout.addWidget(self.create_btn)

        self.check_btn = QPushButton("Check Observable Objects")
        self.check_btn.clicked.connect(self.check_observability)
        left_layout.addWidget(self.check_btn)

        left_layout.addStretch()
        main_layout.addWidget(left_widget, 1)

        # Right Container (Table)
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(6)
        self.results_table.setHorizontalHeaderLabels(["Name", "Type", "RA", "Dec", "Transit", "Mag"])
        
        # Add tooltips to headers to explain columns
        self.results_table.horizontalHeaderItem(0).setToolTip("Primary Identifier")
        self.results_table.horizontalHeaderItem(1).setToolTip("Object Classification")
        self.results_table.horizontalHeaderItem(2).setToolTip("Right Ascension (ICRS)")
        self.results_table.horizontalHeaderItem(3).setToolTip("Declination (ICRS)")
        self.results_table.horizontalHeaderItem(4).setToolTip("Time of Meridian Transit (Highest point in the sky)")
        self.results_table.horizontalHeaderItem(5).setToolTip("Visual Magnitude (Lower is brighter)")
        
        main_layout.addWidget(self.results_table, 2)

        # Load initial config
        self.load_config()
        self.load_session()

    def load_config(self):
        """
        Reads 'config.csv' to pre-fill the observer's location (Lat, Lon, Elev, Name).
        """
        try:
            with open('config.csv', newline='') as csvfile:
                c = csv.reader(csvfile)
                next(c)  # Skip header
                row = next(c)
                if len(row) == 1:
                    row = row[0].split(',')
                
                lat, lon, elev, name = row
                self.lat_edit.setText(lat.strip())
                self.lon_edit.setText(lon.strip())
                self.elev_edit.setText(elev.strip())
                self.name_edit.setText(name.strip())
                self.log("Configuration loaded from config.csv")
        except Exception as e:
            self.log(f"Config load info: {e}")

    def load_session(self):
        """Loads the last used start and end times from session.json."""
        try:
            if os.path.exists('session.json'):
                with open('session.json', 'r') as f:
                    data = json.load(f)
                    if 'start' in data: self.time_begin_edit.setText(data['start'])
                    if 'end' in data: self.time_end_edit.setText(data['end'])
                    self.log("Restored last session times.")
        except Exception as e:
            self.log(f"Could not load session: {e}")

    def save_session(self):
        """Saves the current start and end times to session.json."""
        try:
            data = {
                'start': self.time_begin_edit.text(),
                'end': self.time_end_edit.text()
            }
            with open('session.json', 'w') as f:
                json.dump(data, f)
        except Exception as e:
            self.log(f"Could not save session: {e}")

    def log(self, message):
        """Helper to print to console and update the GUI status bar."""
        print(message)
        self.statusBar().showMessage(message)

    def create_observer(self):
        """
        Instantiates the astroplan.Observer object based on user input.
        This object is required for all subsequent astronomical calculations.
        """
        try:
            lat = float(self.lat_edit.text()) * u.deg
            lon = float(self.lon_edit.text()) * u.deg
            elev = float(self.elev_edit.text()) * u.m
            name = self.name_edit.text()

            tf = TimezoneFinder()
            # Automatically determine timezone string (e.g., 'America/New_York') from coordinates
            timezone_str = tf.timezone_at(lng=lon.value, lat=lat.value)
            
            # Create the Observer object
            self.observer = astroplan.Observer(latitude=lat, longitude=lon, elevation=elev, name=name, timezone=timezone_str)

            if self.temp_edit.text(): self.observer.temperature = float(self.temp_edit.text()) * u.deg_C
            if self.press_edit.text(): self.observer.pressure = float(self.press_edit.text()) * u.mbar
            if self.hum_edit.text(): self.observer.relative_humidity = float(self.hum_edit.text()) / 100.0

            self.log(f"\nObserver Created: {self.observer.name}\nLocation: {self.observer.location}\nTimezone: {self.observer.timezone}")
            if self.time_begin_edit.text(): self.log(f"Start: {self.time_begin_edit.text()}")
            if self.time_end_edit.text(): self.log(f"End: {self.time_end_edit.text()}")

        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))

    def check_observability(self):
        """
        Performs a grid search of the visible sky to find observable objects.
        1. Generates a grid of Alt/Az points representing the local sky.
        2. Queries SIMBAD for objects around those points that meet magnitude criteria.
        """
        if not hasattr(self, 'observer'):
            QMessageBox.warning(self, "Warning", "Please create an observer first.")
            return

        # Get Time
        start_str = self.time_begin_edit.text()
        end_str = self.time_end_edit.text()

        if not start_str or not end_str:
            QMessageBox.warning(self, "Warning", "Please enter start and end times.")
            return

        try:
            t_start = Time(start_str)
            
            # Save times for next session
            self.save_session()
            
            # Get Horizon Limit
            horizon = 0.0
            if self.horizon_edit.text():
                horizon = float(self.horizon_edit.text())

            # Get Magnitude Limit
            mag_limit = 6.0
            if self.mag_edit.text():
                mag_limit = float(self.mag_edit.text())

            # Configure SIMBAD
            custom_simbad = Simbad()
            custom_simbad.ROW_LIMIT = 0

            # Define Grid of Cone Queries
            # Strategy: Instead of querying the whole sky, we query specific regions 
            # that are currently visible to the observer.
            # We start 10 degrees above the user-defined horizon limit to ensure good visibility.
            alt_min = horizon + 10.0
            if alt_min >= 90: alt_min = 89.0
            
            # Generate grid points (Alt/Az)
            grid_coords = []
            step = 15 # Step size in degrees for the grid
            for alt in range(int(alt_min), 90, step):
                for az in range(0, 360, step):
                    # Convert local AltAz coordinates to Celestial ICRS (RA/Dec) at the specific start time
                    c = SkyCoord(alt=alt*u.deg, az=az*u.deg, frame='altaz', 
                                 obstime=t_start, location=self.observer.location)
                    grid_coords.append(c.transform_to('icrs'))

            self.log(f"Starting SIMBAD grid search ({len(grid_coords)} queries)...")
            self.log("Scanning region >10 deg above horizon. Please wait...")
            
            unique_objects = {} # Dictionary to store results and prevent duplicates (Key: Main ID)
            
            for i, coord in enumerate(grid_coords):
                self.log(f"Querying region {i+1}/{len(grid_coords)}...")
                # Keep GUI responsive
                QApplication.processEvents()
                
                try:
                    # Server-side filtering using TAP (ADQL)
                    # We use ADQL (Astronomical Data Query Language) to filter data on the server.
                    # This is much faster than downloading everything and filtering in Python.
                    ra_deg = coord.ra.deg
                    dec_deg = coord.dec.deg
                    
                    # Query Explanation:
                    # JOIN basic and flux tables.
                    # CONTAINS/POINT/CIRCLE: Select objects within 10 degrees of our grid point.
                    # flux."filter" = 'V': Look for Visual magnitude data.
                    # flux.flux < {mag_limit}: Filter for brightness (lower magnitude = brighter).
                    query = f"""
                        SELECT basic.main_id, basic.ra, basic.dec, basic.otype, flux.flux
                        FROM basic
                        JOIN flux ON basic.oid = flux.oidref
                        WHERE 1=CONTAINS(POINT('ICRS', basic.ra, basic.dec), CIRCLE('ICRS', {ra_deg}, {dec_deg}, 10))
                        AND flux."filter" = 'V'
                        AND flux.flux < {mag_limit}
                    """
                    
                    result_table = custom_simbad.query_tap(query)
                    
                    if result_table:
                        for row in result_table:
                            try:
                                # Decode bytes to string if necessary (common issue with VOTable data)
                                name = row['main_id']
                                if isinstance(name, bytes): name = name.decode('utf-8')
                                name = str(name)
                                
                                # Skip if we already found this object in a previous grid overlap
                                if name in unique_objects:
                                    continue
                                
                                ra = float(row['ra'])
                                dec = float(row['dec'])
                                otype = row['otype']
                                if isinstance(otype, bytes): otype = otype.decode('utf-8')
                                otype = str(otype).strip()
                                otype = OTYPE_MAP.get(otype, otype) # Convert to human readable
                                mag = float(row['flux'])
                                
                                unique_objects[name] = {
                                    'name': name,
                                    'type': otype,
                                    'ra': ra,
                                    'dec': dec,
                                    'mag': mag,
                                    'target': FixedTarget(coord=SkyCoord(ra=ra*u.deg, dec=dec*u.deg), name=name)
                                }
                            except Exception:
                                continue
                except Exception as e:
                    print(f"Query failed: {e}", file=sys.stderr)
                    self.log(f"Query failed: {e}")
                
                # 1 second delay to respect SIMBAD server usage policies
                time.sleep(1)

            # Display Results
            self.results_table.setSortingEnabled(False)
            self.results_table.setRowCount(0)

            for obj in unique_objects.values():
                # Calculate the time the object crosses the meridian (highest point in sky)
                try:
                    transit = self.observer.target_meridian_transit_time(t_start, obj['target'], which='next')
                    transit_str = transit.iso.split(' ')[1][:5]
                except:
                    transit_str = "-"

                row = self.results_table.rowCount()
                self.results_table.insertRow(row)
                self.results_table.setItem(row, 0, QTableWidgetItem(obj['name']))
                self.results_table.setItem(row, 1, QTableWidgetItem(obj['type']))
                self.results_table.setItem(row, 2, NumericTableWidgetItem(f"{obj['ra']:.4f}"))
                self.results_table.setItem(row, 3, NumericTableWidgetItem(f"{obj['dec']:.4f}"))
                self.results_table.setItem(row, 4, QTableWidgetItem(transit_str))
                self.results_table.setItem(row, 5, NumericTableWidgetItem(f"{obj['mag']:.2f}"))
            
            self.results_table.setSortingEnabled(True)
            self.log(f"Found {len(unique_objects)} visible objects.")

        except Exception as e:
            self.log(f"Error checking observability: {e}")

if __name__ == "__main__":
    # Download IERS data (Earth rotation data) required for precise time/coordinate conversions
    print("Checking IERS data...")
    try:
        astroplan.download_IERS_A(show_progress=False)
    except:
        pass
    app = QApplication(sys.argv)
    window = AdAstraWindow()
    window.show()
    sys.exit(app.exec())