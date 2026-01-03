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
from datetime import datetime
import pytz
import numpy as np
import astropy.units as u
import astroplan
from astroplan import FixedTarget, AltitudeConstraint, AtNightConstraint
from astropy.time import Time
from timezonefinder import TimezoneFinder
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from PyQt6.QtGui import QIcon
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                             QLineEdit, QPushButton, QTextEdit, QFormLayout, QMessageBox,
                             QTableWidget, QTableWidgetItem, QComboBox, QLabel, QDialog, 
                             QListWidget)

# Mapping of SIMBAD O-types to human-readable descriptions
# Source: https://simbad.cds.unistra.fr/guide/otypes.htx
OTYPE_MAP = {
     '?': 'Object of Unknown Nature',
     '..1': '{pr*} Pre-Main Sequence Star ',
     '..10': 'Barium Star',
     '..11': 'Dwarf Carbon Star',
     '..12': 'Carbon-Enhanced Metal Poor Star',
     '..13': '{Al*} Eclipsing Binary of Algol type',
     '..14': '{bL*}Eclipsing Binary of beta Lyr type',
     '..15': '{WU*} Eclipsing Binary of W UMa type',
     '..16': '{NL*} Nova-like Binary',
     '..17': '{DN*} Dwarf Nova',
     '..18': '{DQ*} CV of DQ Her type  Intermediate polar.',
     '..19': '{AM*} CV of AM CVn type',
     '..2': 'LBV=Luminous Blue Variable',
     '..20': 'Irregular Variable with rapid variations',
     '..21': '{Fl*} Flare Star',
     '..22': 'Star showing Eclipses by its Planet',
     '..23': '{*iC} Star towards a Cluster',
     '..24': '{*iA} Star towards an Association',
     '..25': '{*iN} Star towards a Nebula',
     '..26': '{*i*} Star in double system',
     '..27': '{BNe} Bright Nebula',
     '..28': '{HzG} Galaxy with high redshift',
     '..29': '{ERO} ERO/VRO, Extremely/Very Red Object',
     '..3': '{FU*} FU Ori Variable',
     '..30': 'ULIRG, Ultra Luminous Infrared Galaxy',
     '..31': '{LyA, DLA, mAL, LLS, BAL} Absorption Line System',
     '..32': '{red} Very Red Source',
     '..4': 'Red Clump Star',
     '..5': '{sr*} Semi-Regular Variable',
     '..6': 'O-rich AGB Star',
     '..7': '{ZZ*} Pulsating White Dwarf',
     '..8': 'ELMWD=Extremely Low Mass White Dwarf',
     '..9': 'CH Star',
     '*': 'Star',
     '**': 'Double or Multiple Star',
     'a2*': 'alpha2 CVn Variable',
     'AB*': 'Asymptotic Giant Branch Star',
     'Ae*': 'Herbig Ae/Be Star',
     'AGN': 'Active Galaxy Nucleus',
     'As*': 'Association of Stars',
     'bC*': 'beta Cep Variable',
     'bCG': 'Blue Compact Galaxy',
     'BD*': 'Brown Dwarf',
     'Be*': 'Be Star',
     'BH': 'Black Hole',
     'BiC': 'Brightest Galaxy in a Cluster (BCG)',
     'Bla': 'Blazar',
     'BLL': 'BL Lac',
     'blu': 'Blue Object',
     'BS*': 'Blue Straggler',
     'bub': 'Bubble',
     'BY*': 'BY Dra Variable',
     'C*': 'Carbon Star',
     'cC*': 'Classical Cepheid Variable',
     'Ce*': 'Cepheid Variable',
     'CGb': 'Cometary Globule / Pillar',
     'CGG': 'Compact Group of Galaxies',
     'Cl*': 'Cluster of Stars',
     'Cld': 'Cloud',
     'ClG': 'Cluster of Galaxies',
     'cm': 'Centimetric Radio Source',
     'cor': 'Dense Core',
     'CV*': 'Cataclysmic Binary',
     'DNe': 'Dark Cloud (nebula)',
     'dS*': 'delta Sct Variable',
     'EB*': 'Eclipsing Binary',
     'El*': 'Ellipsoidal Variable',
     'Em*': 'Emission-line Star',
     'EmG': 'Emission-line galaxy',
     'EmO': 'Emission Object',
     'Er*': 'Eruptive Variable',
     'err': 'Not an Object (Error, Artefact, ...)',
     'ev': 'Transient Event',
     'Ev*': 'Evolved Star',
     'FIR': 'Far-IR source (λ >= 30 µm)',
     'flt': 'Interstellar Filament',
     'G': 'Galaxy',
     'gam': 'Gamma-ray Source',
     'gB': 'Gamma-ray Burst',
     'gD*': 'gamma Dor Variable',
     'GiC': 'Galaxy towards a Cluster of Galaxies',
     'GiG': 'Galaxy towards a Group of Galaxies',
     'GiP': 'Galaxy in Pair of Galaxies',
     'glb': 'Globule (low-mass dark cloud)',
     'GlC': 'Globular Cluster',
     'gLe': 'Gravitational Lens',
     'gLS': 'Gravitational Lens System (lens+images)',
     'GNe': 'Nebula',
     'GrG': 'Group of Galaxies',
     'grv': 'Gravitational Source',
     'GWE': 'Gravitational Wave Event',
     'H2G': 'HII Galaxy',
     'HB*': 'Horizontal Branch Star',
     'HH': 'Herbig-Haro Object',
     'HI': 'HI (21cm) Source',
     'HII': 'HII Region',
     'HS*': 'Hot Subdwarf',
     'HV*': 'High Velocity Star',
     'HVC': 'High-velocity Cloud',
     'HXB': 'High Mass X-ray Binary',
     'IG': 'Interacting Galaxies',
     'IR': 'Infra-Red Source',
     'Ir*': 'Irregular Variable',
     'ISM': 'Interstellar Medium Object',
     'LeG': 'Gravitationally Lensed Image of a Galaxy',
     'LeI': 'Gravitationally Lensed Image',
     'LeQ': 'Gravitationally Lensed Image of a Quasar',
     'Lev': '(Micro)Lensing Event',
     'LIN': 'LINER-type Active Galaxy Nucleus',
     'LM*': 'Low-mass Star',
     'LP*': 'Long-Period Variable',
     'LSB': 'Low Surface Brightness Galaxy',
     'LXB': 'Low Mass X-ray Binary',
     'Ma*': 'Massive Star',
     'Mas': 'Maser',
     'MGr': 'Moving Group',
     'Mi*': 'Mira Variable',
     'MIR': 'Mid-IR Source (3 to 30 µm)',
     'mm': 'Millimetric Radio Source',
     'MoC': 'Molecular Cloud',
     'mR': 'Metric Radio Source',
     'MS*': 'Main Sequence Star',
     'mul': 'Composite Object, Blend',
     'N*': 'Neutron Star',
     'NIR': 'Near-IR Source (λ < 3 µm)',
     'No*': 'Classical Nova',
     'OH*': 'OH/IR Star',
     'OpC': 'Open Cluster',
     'Opt': 'Optical Source',
     'Or*': 'Orion Variable',
     'out': 'Outflow',
     'pA*': 'Post-AGB Star',
     'PaG': 'Pair of Galaxies',
     'PCG': 'Proto Cluster of Galaxies',
     'Pe*': 'Chemically Peculiar Star',
     'Pl': 'Extra-solar Planet',
     'PM*': 'High Proper Motion Star',
     'PN': 'Planetary Nebula',
     'PoC': 'Part of Cloud',
     'PoG': 'Part of a Galaxy',
     'Psr': 'Pulsar',
     'Pu*': 'Pulsating Variable',
     'QSO': 'Quasar',
     'Rad': 'Radio Source',
     'rB': 'Radio Burst',
     'RC*': 'R CrB Variable',
     'reg': 'Region defined in the Sky',
     'rG': 'Radio Galaxy',
     'RG*': 'Red Giant Branch star',
     'RNe': 'Reflection Nebula',
     'Ro*': 'Rotating Variable',
     'RR*': 'RR Lyrae Variable',
     'RS*': 'RS CVn Variable',
     'RV*': 'RV Tauri Variable',
     'S*': 'S Star',
     's*b': 'Blue Supergiant',
     's*r': 'Red Supergiant',
     's*y': 'Yellow Supergiant',
     'SB*': 'Spectroscopic Binary',
     'SBG': 'Starburst Galaxy',
     'SCG': 'Supercluster of Galaxies',
     'SFR': 'Star Forming Region',
     'sg*': 'Evolved Supergiant',
     'sh': 'Interstellar Shell',
     'smm': 'Sub-Millimetric Source',
     'SN*': 'SuperNova',
     'SNR': 'SuperNova Remnant',
     'St*': 'Stellar Stream',
     'SX*': 'SX Phe Variable',
     'Sy*': 'Symbiotic Star',
     'Sy1': 'Seyfert 1 Galaxy',
     'Sy2': 'Seyfert 2 Galaxy',
     'SyG': 'Seyfert Galaxy',
     'TT*': 'T Tauri Star',
     'ULX': 'Ultra-luminous X-ray Source',
     'UV': 'UV-emission Source',
     'V*': 'Variable Star',
     'var': 'Variable source',
     'vid': 'Underdense Region of the Universe',
     'WD*': 'White Dwarf',
     'WR*': 'Wolf-Rayet',
     'WV*': 'Type II Cepheid Variable',
     'X': 'X-ray Source',
     'XB*': 'X-ray Binary',
     'Y*O': 'Young Stellar Object',
 }

class NumericTableWidgetItem(QTableWidgetItem):
    """
    Custom TableItem to ensure numerical sorting for columns like Magnitude and Coordinates.
    Standard QTableWidgetItem sorts alphabetically (e.g., "10" comes before "2").
    """
    def __init__(self, text, sort_value=None):
        super().__init__(text)
        self.sort_value = sort_value

    def __lt__(self, other):
        if (getattr(self, 'sort_value', None) is not None and 
            getattr(other, 'sort_value', None) is not None):
            return self.sort_value < other.sort_value

        try:
            return float(self.text()) < float(other.text())
        except ValueError:
            return super().__lt__(other)

class LocationManagerDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Manage Locations")
        self.resize(500, 400)
        self.parent_window = parent
        
        layout = QVBoxLayout(self)
        
        # List of locations
        layout.addWidget(QLabel("Saved Locations:"))
        self.loc_list = QListWidget()
        self.loc_list.addItems(sorted(self.parent_window.locations.keys()))
        self.loc_list.itemClicked.connect(self.load_location)
        layout.addWidget(self.loc_list)
        
        # Form
        form_layout = QFormLayout()
        self.name_edit = QLineEdit()
        self.lat_edit = QLineEdit()
        self.lon_edit = QLineEdit()
        self.elev_edit = QLineEdit()
        
        form_layout.addRow("Name:", self.name_edit)
        form_layout.addRow("Latitude (deg):", self.lat_edit)
        form_layout.addRow("Longitude (deg):", self.lon_edit)
        form_layout.addRow("Elevation (m):", self.elev_edit)
        layout.addLayout(form_layout)
        
        # Buttons
        btn_layout = QHBoxLayout()
        self.new_btn = QPushButton("New / Clear")
        self.new_btn.clicked.connect(self.clear_form)
        self.save_btn = QPushButton("Save")
        self.save_btn.clicked.connect(self.save_location)
        self.del_btn = QPushButton("Delete")
        self.del_btn.clicked.connect(self.delete_location)
        
        btn_layout.addWidget(self.new_btn)
        btn_layout.addWidget(self.save_btn)
        btn_layout.addWidget(self.del_btn)
        layout.addLayout(btn_layout)

    def load_location(self, item):
        name = item.text()
        if name in self.parent_window.locations:
            data = self.parent_window.locations[name]
            self.name_edit.setText(name)
            self.lat_edit.setText(str(data['latitude']))
            self.lon_edit.setText(str(data['longitude']))
            self.elev_edit.setText(str(data['elevation']))

    def clear_form(self):
        self.loc_list.clearSelection()
        self.name_edit.clear()
        self.lat_edit.clear()
        self.lon_edit.clear()
        self.elev_edit.clear()

    def save_location(self):
        name = self.name_edit.text().strip()
        if not name:
            QMessageBox.warning(self, "Error", "Name cannot be empty.")
            return
        try:
            lat = float(self.lat_edit.text())
            lon = float(self.lon_edit.text())
            elev = float(self.elev_edit.text())
            
            self.parent_window.locations[name] = {
                'latitude': lat, 'longitude': lon, 'elevation': elev
            }
            self.parent_window.save_data()
            self.refresh_list()
            self.parent_window.refresh_locations()
            QMessageBox.information(self, "Success", f"Location '{name}' saved.")
        except ValueError:
            QMessageBox.warning(self, "Error", "Latitude, Longitude, and Elevation must be valid numbers.")

    def delete_location(self):
        name = self.name_edit.text().strip()
        if not name: return
        
        if name in self.parent_window.locations:
            reply = QMessageBox.question(self, 'Confirm Delete', f"Delete '{name}'?",
                                         QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            if reply == QMessageBox.StandardButton.Yes:
                del self.parent_window.locations[name]
                self.parent_window.save_data()
                self.clear_form()
                self.refresh_list()
                self.parent_window.refresh_locations()
        else:
            QMessageBox.warning(self, "Error", "Location not found.")

    def refresh_list(self):
        self.loc_list.clear()
        self.loc_list.addItems(sorted(self.parent_window.locations.keys()))

class AdAstraWindow(QMainWindow):
    def __init__(self):
        """
        Initialize the Main Window, Layouts, and Widgets.
        """
        super().__init__()
        
        # In-memory store for location data
        self.locations = {}
        
        self.setWindowTitle("Ad Astra")
        self.resize(1000, 600)
        
        # Set window icon (requires 'icon.png' in the same directory)
        self.setWindowIcon(QIcon('icon.png'))

        # Main widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)

        # Left Container
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)

        # Form Layout for inputs
        form_layout = QFormLayout()

        # Location selection
        self.location_combo = QComboBox()
        self.location_combo.currentIndexChanged.connect(self.on_location_change)
        
        self.new_loc_btn = QPushButton("New location")
        self.new_loc_btn.clicked.connect(self.open_location_manager)

        loc_layout = QHBoxLayout()
        loc_layout.addWidget(self.location_combo)
        loc_layout.addWidget(self.new_loc_btn)
        
        form_layout.addRow("Location:", loc_layout)
        
        # Plain text details
        self.lat_label = QLabel("-")
        self.lon_label = QLabel("-")
        self.elev_label = QLabel("-")
        
        form_layout.addRow("Latitude:", self.lat_label)
        form_layout.addRow("Longitude:", self.lon_label)
        form_layout.addRow("Elevation:", self.elev_label)

        # User inputs
        self.time_begin_edit = QLineEdit()
        self.time_begin_edit.setPlaceholderText("YYYY-MM-DD HH:MM")
        form_layout.addRow("Start Time:", self.time_begin_edit)

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
        self.create_btn = QPushButton("Create Location")
        self.create_btn.clicked.connect(self.toggle_observer_state)
        left_layout.addWidget(self.create_btn)

        self.check_btn = QPushButton("Check Observable Objects")
        self.check_btn.clicked.connect(self.check_observability)
        self.check_btn.setVisible(False)
        left_layout.addWidget(self.check_btn)

        left_layout.addStretch()
        main_layout.addWidget(left_widget, 1)

        # Right Container (Table)
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(6)
        self.results_table.setHorizontalHeaderLabels(["Name", "Type", "RA", "Dec", "Duration", "Mag"])
        
        # Add tooltips to headers to explain columns
        self.results_table.horizontalHeaderItem(0).setToolTip("Primary Identifier")
        self.results_table.horizontalHeaderItem(1).setToolTip("Object Classification")
        self.results_table.horizontalHeaderItem(2).setToolTip("Right Ascension (ICRS)")
        self.results_table.horizontalHeaderItem(3).setToolTip("Declination (ICRS)")
        self.results_table.horizontalHeaderItem(4).setToolTip("Duration of visibility above horizon during observation window")
        self.results_table.horizontalHeaderItem(5).setToolTip("Visual Magnitude (Lower is brighter)")
        
        main_layout.addWidget(self.results_table, 2)

        # Load initial config
        self.load_data()

    def load_data(self):
        """Loads all locations and session data from adastra_data.json."""
        data_file = 'adastra_data.json'
        if not os.path.exists(data_file):
            self.log("No data file found. Trying to import from old config.csv...")
            self.locations = {}
            try:
                with open('config.csv', newline='') as csvfile:
                    c = csv.reader(csvfile)
                    next(c)  # Skip header
                    row = next(c)
                    if len(row) == 1: row = row[0].split(',')
                    lat, lon, elev, name = row
                    self.locations[name.strip()] = {
                        "latitude": float(lat), "longitude": float(lon), "elevation": float(elev)
                    }
                    self.log("Successfully imported location from config.csv.")
                    self.save_data() # Create the new data file
            except Exception as e:
                self.log(f"Could not import from config.csv: {e}")

        try:
            with open(data_file, 'r') as f:
                data = json.load(f)
                self.locations = data.get('locations', {})
                
                # Populate dropdown
                self.location_combo.blockSignals(True)
                self.location_combo.clear()
                self.location_combo.addItems(sorted(self.locations.keys()))
                self.location_combo.blockSignals(False)

                # Set to last used location
                last_loc = data.get('last_used_location')
                if last_loc in self.locations:
                    self.location_combo.setCurrentText(last_loc)
                
                # Trigger update of labels
                self.on_location_change(self.location_combo.currentIndex())
                
                # Restore times
                self.time_begin_edit.setText(data.get('last_start_time', ''))
                self.log("Loaded saved locations and session data.")

        except Exception as e:
            self.log(f"Could not load data file: {e}")

    def save_data(self):
        """Saves all locations and session data to adastra_data.json."""
        try:
            data = {
                'last_used_location': self.location_combo.currentText(),
                'last_start_time': self.time_begin_edit.text(),
                'locations': self.locations
            }
            with open('adastra_data.json', 'w') as f:
                json.dump(data, f, indent=4)
            self.log("Saved locations and session data.")
        except Exception as e:
            self.log(f"Could not save data: {e}")

    def on_location_change(self, index):
        """Updates the Lat/Lon/Elev fields when a new location is selected."""
        loc_name = self.location_combo.currentText()
        if loc_name in self.locations:
            loc_data = self.locations[loc_name]
            self.lat_label.setText(str(loc_data.get('latitude', '-')))
            self.lon_label.setText(str(loc_data.get('longitude', '-')))
            self.elev_label.setText(str(loc_data.get('elevation', '-')))
        else:
            self.lat_label.setText("-")
            self.lon_label.setText("-")
            self.elev_label.setText("-")

    def open_location_manager(self):
        dialog = LocationManagerDialog(self)
        dialog.exec()

    def refresh_locations(self):
        current = self.location_combo.currentText()
        self.location_combo.blockSignals(True)
        self.location_combo.clear()
        self.location_combo.addItems(sorted(self.locations.keys()))
        if current in self.locations:
            self.location_combo.setCurrentText(current)
        self.location_combo.blockSignals(False)
        self.on_location_change(self.location_combo.currentIndex())

    def log(self, message):
        """Helper to print to console and update the GUI status bar."""
        print(message)
        self.statusBar().showMessage(message)

    def toggle_observer_state(self):
        """Toggles between Setup mode (inputs enabled) and Ready mode (inputs disabled)."""
        if self.check_btn.isVisible():
            # Currently in Ready mode -> Switch to Setup
            self.set_inputs_enabled(True)
            self.check_btn.setVisible(False)
            self.create_btn.setText("Update Observer")
        else:
            # Currently in Setup mode -> Switch to Ready
            if self.create_observer():
                self.set_inputs_enabled(False)
                self.check_btn.setVisible(True)
                self.create_btn.setText("Change Settings")

    def set_inputs_enabled(self, enabled):
        """Enables or disables configuration inputs."""
        widgets = [
            self.location_combo, self.new_loc_btn,
            self.time_begin_edit,
            self.horizon_edit, self.mag_edit,
            self.temp_edit, self.press_edit, self.hum_edit
        ]
        for w in widgets:
            w.setEnabled(enabled)

    def create_observer(self):
        """
        Instantiates the astroplan.Observer object based on user input.
        This object is required for all subsequent astronomical calculations.
        Returns True if successful, False otherwise.
        """
        try:
            loc_name = self.location_combo.currentText()
            if not loc_name or loc_name not in self.locations:
                raise ValueError("Please select a valid location.")
            loc_data = self.locations[loc_name]
            lat = float(loc_data['latitude']) * u.deg
            lon = float(loc_data['longitude']) * u.deg
            elev = float(loc_data['elevation']) * u.m
            name = loc_name

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
            return True

        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))
            return False

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

        if not start_str:
            QMessageBox.warning(self, "Warning", "Please enter start time.")
            return

        try:
            # Parse start time as local time using the observer's timezone
            dt_naive = datetime.strptime(start_str, "%Y-%m-%d %H:%M")
            tz = self.observer.timezone
            if isinstance(tz, str):
                tz = pytz.timezone(tz)
            dt_aware = tz.localize(dt_naive)
            t_start = Time(dt_aware)
            
            # Automatically calculate end time as next nautical dawn
            t_end = self.observer.twilight_morning_nautical(t_start, which='next')
            self.log(f"Observation window: {t_start.iso} to {t_end.iso} (Nautical Dawn)")

            # Get Horizon Limit
            horizon = 0.0
            if self.horizon_edit.text():
                horizon = float(self.horizon_edit.text())

            # Get Magnitude Limit
            mag_limit = 6.0
            if self.mag_edit.text():
                mag_limit = float(self.mag_edit.text())
            
            # Save all current settings for the next session
            self.save_data()

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

            # Calculate visibility duration for all objects
            objects_list = list(unique_objects.values())
            targets = [obj['target'] for obj in objects_list]
            duration_data = []

            if targets and t_end > t_start:
                try:
                    time_step_mins = 5
                    time_step = time_step_mins * u.min
                    delta_time = (t_end - t_start).to(u.min).value
                    n_steps = int(delta_time / time_step_mins) + 1
                    time_grid = t_start + (np.arange(n_steps) * time_step)
                    
                    altaz_grid = self.observer.altaz(time_grid, targets, grid_times_targets=True)
                    is_above_horizon = altaz_grid.alt > (horizon * u.deg)
                    
                    # Since t_end is fixed to Nautical Dawn, we don't need to filter for night again.
                    is_visible = is_above_horizon
                    
                    visible_counts = np.sum(is_visible, axis=1)
                    
                    for count in visible_counts:
                        mins_total = int(count * time_step_mins)
                        h = mins_total // 60
                        m = mins_total % 60
                        duration_data.append((f"{h}h {m}m", mins_total))
                except Exception as e:
                    print(f"Duration calculation error: {e}", file=sys.stderr)
                    self.log(f"Duration calculation failed: {e}")
                    duration_data = [("-", 0)] * len(targets)
            else:
                duration_data = [("-", 0)] * len(targets)

            for i, obj in enumerate(objects_list):
                duration_str, duration_mins = duration_data[i]
                
                row = self.results_table.rowCount()
                self.results_table.insertRow(row)
                self.results_table.setItem(row, 0, QTableWidgetItem(obj['name']))
                self.results_table.setItem(row, 1, QTableWidgetItem(obj['type']))
                
                ra_str = obj['target'].coord.ra.to_string(unit=u.hour, sep=('h', 'm', 's'), precision=1, pad=True)
                dec_str = obj['target'].coord.dec.to_string(unit=u.deg, sep=('d', 'm', 's'), precision=1, alwayssign=True, pad=True)
                
                self.results_table.setItem(row, 2, NumericTableWidgetItem(ra_str, sort_value=obj['ra']))
                self.results_table.setItem(row, 3, NumericTableWidgetItem(dec_str, sort_value=obj['dec']))
                self.results_table.setItem(row, 4, NumericTableWidgetItem(duration_str, sort_value=duration_mins))
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
    app.setWindowIcon(QIcon('icon.png'))
    window = AdAstraWindow()
    window.show()
    sys.exit(app.exec())