from flask import Flask, render_template, request, jsonify, redirect, url_for
import io
import base64
from matplotlib.figure import Figure
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astroplan import Observer, FixedTarget, AtNightConstraint, AltitudeConstraint
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astroquery.hips2fits import hips2fits
from astroquery.ipac.ned import Ned
from timezonefinder import TimezoneFinder
from datetime import datetime
import pytz

# Import our new Gaia module
from gaia_module import query_gaia_cone

app = Flask(__name__)
app.title = "AstroSlurp"

# --- Configuration ---
PLOT_TEXT_COLOR = 'white'
PLOT_BG_COLOR = '#2b2b35'

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

# --- Helper Functions ---

def create_observer(lat, lon, elev):
    """Creates an astroplan Observer object."""
    tf = TimezoneFinder()
    timezone_str = tf.timezone_at(lng=lon, lat=lat)
    observer = Observer(latitude=lat*u.deg, longitude=lon*u.deg, elevation=elev*u.m, timezone=timezone_str)
    return observer

def get_deep_dive_data(object_name):
    """Fetches detailed data for a single object from SIMBAD."""
    try:
        # Request multiple flux fields to ensure we get a magnitude
        simbad = Simbad()
        simbad.add_votable_fields('plx', 'otype', 'flux(V)', 'flux(G)', 'flux(R)', 'flux(B)', 'rv_value')
        table = simbad.query_object(object_name)
        
    except Exception as e:
        # Fallback if specific fields cause error (though add_votable_fields usually handles missing ones gracefully)
        print(f"Primary SIMBAD query failed: {e}. Retrying with minimal fields.")
        try:
            simbad = Simbad()
            simbad.add_votable_fields('plx', 'otype')
            table = simbad.query_object(object_name)
        except Exception as e2:
             return None, f"Secondary SIMBAD query failed: {e2}"

    if not table:
        return None, "Object not found in SIMBAD."

    row = table[0]
    
    # Robust Magnitude Logic (Waterfall)
    mag = "-"
    mag_type = ""
    
    # Check V (Visual)
    if 'FLUX_V' in row.colnames and not np.ma.is_masked(row['FLUX_V']):
        mag = f"{row['FLUX_V']:.2f}"
        mag_type = "(V)"
    # Check G (Gaia) - often available for stars where V is missing
    elif 'FLUX_G' in row.colnames and not np.ma.is_masked(row['FLUX_G']):
        mag = f"{row['FLUX_G']:.2f}"
        mag_type = "(G)"
    # Check R (Red)
    elif 'FLUX_R' in row.colnames and not np.ma.is_masked(row['FLUX_R']):
        mag = f"{row['FLUX_R']:.2f}"
        mag_type = "(R)"
    # Check B (Blue)
    elif 'FLUX_B' in row.colnames and not np.ma.is_masked(row['FLUX_B']):
        mag = f"{row['FLUX_B']:.2f}"
        mag_type = "(B)"

    # Get object type and convert to human-readable name
    otype_code = row['OTYPE']
    otype_readable = OTYPE_MAP.get(otype_code, otype_code)

    data = {
        'name': row['MAIN_ID'],
        'type': otype_readable,
        'ra': row['RA'],
        'dec': row['DEC'],
        'mag': f"{mag} {mag_type}",
        'coord': SkyCoord(ra=row['RA'], dec=row['DEC'], unit=(u.hourangle, u.deg), frame='icrs')
    }
    
    # Calculate distance
    distance_str = '-'
    
    # 1. Try SIMBAD Parallax
    if 'PLX_VALUE' in row.colnames and not np.ma.is_masked(row['PLX_VALUE']) and row['PLX_VALUE'] > 0:
        d_pc = 1000.0 / row['PLX_VALUE']
        d_ly = d_pc * 3.26156
        distance_str = f"{d_ly:,.1f} ly"
    
    # 2. Fallback to Gaia if SIMBAD failed
    if distance_str == '-':
        try:
            # Query Gaia using the coordinates from SIMBAD
            # Use a small radius (e.g., 2 arcseconds) to find the exact match
            gaia_results = query_gaia_cone(data['coord'].ra.deg, data['coord'].dec.deg, radius_deg=0.0005, limit=1)
            
            if gaia_results and gaia_results[0]['parallax'] is not None and gaia_results[0]['parallax'] > 0:
                plx = gaia_results[0]['parallax']
                d_pc = 1000.0 / plx
                d_ly = d_pc * 3.26156
                distance_str = f"{d_ly:,.1f} ly (Gaia)"
        except Exception as e:
            print(f"Gaia fallback failed: {e}")

    # 3. Fallback to NED (Redshift-Independent Distances) for Galaxies
    if distance_str == '-':
        try:
            # Query NED for redshift-independent distances
            ned_table = Ned.get_table(object_name, table='redshift_independent_distances')
            if ned_table and len(ned_table) > 0:
                # Look for the 'Mean Distance (Mpc)' column. 
                # The column name is often 'Distance (Mpc)' or similar.
                # We will iterate to find a likely candidate.
                dist_col = None
                for col in ned_table.colnames:
                    if 'Distance' in col or 'D (Mpc)' in col:
                        dist_col = col
                        break
                
                if dist_col:
                    # Take the mean of the measurements
                    d_mpc = np.nanmean(ned_table[dist_col])
                    d_ly = d_mpc * 3.26156 * 1e6
                    distance_str = f"{d_ly:,.0f} ly (NED)"
        except Exception as e:
            print(f"NED fallback failed: {e}")

    # 4. Fallback to Redshift (Hubble Law) if everything else failed
    if distance_str == '-':
        if 'RV_VALUE' in row.colnames and not np.ma.is_masked(row['RV_VALUE']):
            v_r = float(row['RV_VALUE']) # Radial velocity in km/s
            # Only apply Hubble law for receding objects (v > 0)
            # and generally only useful for v > 1000 km/s to avoid local peculiar motions
            if v_r > 500: 
                H0 = 70.0 # Hubble constant in km/s/Mpc
                d_mpc = v_r / H0
                d_ly = d_mpc * 3.26156 * 1e6
                distance_str = f"{d_ly:,.0f} ly (Redshift)"

    data['distance'] = distance_str
        
    # Fetch Aliases and find Common Name
    try:
        ids_table = simbad.query_objectids(data['name'])
        aliases = []
        common_name = None
        
        if ids_table:
            for id_row in ids_table:
                val = id_row[0]
                if isinstance(val, bytes): val = val.decode('utf-8')
                alias = ' '.join(str(val).split())
                aliases.append(alias)
                
                if alias.startswith("NAME ") and not common_name:
                    common_name = alias.replace("NAME ", "")
                    
        data['aliases'] = aliases
        data['common_name'] = common_name
        
    except Exception as e:
        print(f"Warning: Could not fetch aliases: {e}")
        data['aliases'] = []
        data['common_name'] = None
            
    return data, None

def generate_image_plot(image_data, title, cmap='gray'):
    """Generates a base64-encoded plot for an image."""
    fig = Figure(figsize=(6, 6), facecolor=PLOT_BG_COLOR)
    ax = fig.add_subplot(111)
    
    h, w = 0, 0
    if image_data.ndim == 3:
        if image_data.shape[0] == 3:
            image_data = np.transpose(image_data, (1, 2, 0))
        norm_data = image_data.astype(float)
        dmin, dmax = np.nanmin(norm_data), np.nanmax(norm_data)
        if dmax > dmin:
            norm_data = (norm_data - dmin) / (dmax - dmin)
        ax.imshow(norm_data, origin='lower')
        h, w, _ = norm_data.shape
    else:
        from astropy.visualization import AsinhStretch, ImageNormalize, simple_norm
        if cmap == 'gray':
            norm = ImageNormalize(stretch=AsinhStretch(a=0.1), vmin=np.nanmin(image_data), vmax=np.nanmax(image_data))
        else: # For scientific data, use a simpler norm
            norm = simple_norm(image_data, 'sqrt', percent=99)
        ax.imshow(image_data, origin='lower', cmap=cmap, norm=norm)
        h, w = image_data.shape
        
    ax.set_title(title, color=PLOT_TEXT_COLOR)
    ax.axis('off')
    
    # Add target marker
    if h > 0 and w > 0:
        ax.plot(w / 2 - 0.5, h / 2 - 0.5, 'o', color='white', markersize=12, alpha=0.7)
        ax.plot(w / 2 - 0.5, h / 2 - 0.5, 'ko', markersize=5)

    fig.tight_layout()
    
    buf = io.BytesIO()
    fig.savefig(buf, format="png", facecolor=fig.get_facecolor())
    data = base64.b64encode(buf.getbuffer()).decode("ascii")
    return f"data:image/png;base64,{data}"

def perform_sky_query(lat, lon, elev, start_time_str, mag_limit, horizon_limit):
    """Performs the main logic for the Sky Query page."""
    try:
        print(f"Starting Sky Query: Lat={lat}, Lon={lon}, MagLimit={mag_limit}")
        observer = create_observer(lat, lon, elev)
        
        # Time window
        dt_naive = datetime.fromisoformat(start_time_str)
        
        # Handle timezone correctly
        if isinstance(observer.timezone, str):
            tz = pytz.timezone(observer.timezone)
        else:
            tz = observer.timezone # It's already a pytz object
            
        dt_aware = tz.localize(dt_naive)
        t_start = Time(dt_aware)
        t_end = observer.twilight_morning_nautical(t_start, which='next')
        
        print(f"Observation window: {t_start.iso} to {t_end.iso}")

        unique_objects = {}
        
        # 1. SIMBAD Query for Deep Sky Objects
        print("Querying SIMBAD for bright deep-sky objects...")
        simbad = Simbad()
        simbad.add_votable_fields('otype', 'flux(V)')
        simbad.TIMEOUT = 120 
        simbad.ROW_LIMIT = 500

        otypes_to_query = ['G', 'Neb', 'OpC', 'GlC', 'SN']
        
        for otype in otypes_to_query:
            try:
                crit = f"Vmag < {mag_limit} & otype = '{otype}'"
                result_table = simbad.query_criteria(crit)
                
                if result_table:
                    for row in result_table:
                        name = row['MAIN_ID']
                        if name not in unique_objects:
                            unique_objects[name] = {
                                'name': name,
                                'type': OTYPE_MAP.get(row['OTYPE'], row['OTYPE']),
                                'ra': row['RA'], 'dec': row['DEC'],
                                'mag': f"{row['FLUX_V']:.2f}" if 'FLUX_V' in row.colnames and not np.ma.is_masked(row['FLUX_V']) else '-',
                                'coord': SkyCoord(ra=row['RA'], dec=row['DEC'], unit=(u.hourangle, u.deg))
                            }
            except Exception as e:
                print(f"SIMBAD criteria query failed for {otype}: {e}")

        # 2. Vizier Query for Bright Stars (Hipparcos Catalog)
        print("Querying Vizier (Hipparcos) for bright stars...")
        try:
            v = Vizier(columns=['HIP', 'RAICRS', 'DEICRS', 'Vmag'], row_limit=-1) # -1 means no limit
            result = v.query_constraints(catalog='I/239/hip_main', Vmag=f"<{mag_limit}")
            
            if result and len(result) > 0:
                table = result[0]
                print(f"Found {len(table)} bright stars from Hipparcos.")
                
                for row in table:
                    hip_id = f"HIP {row['HIP']}"
                    if hip_id not in unique_objects:
                        unique_objects[hip_id] = {
                            'name': hip_id,
                            'type': 'Star',
                            'ra': row['RAICRS'], 'dec': row['DEICRS'], # These are in degrees
                            'mag': f"{row['Vmag']:.2f}",
                            'coord': SkyCoord(ra=row['RAICRS']*u.deg, dec=row['DEICRS']*u.deg)
                        }
        except Exception as e:
            print(f"Vizier query failed: {e}")

        # Calculate observability
        results = []
        print(f"Calculating observability for {len(unique_objects)} candidates...")
        
        frame_start = AltAz(obstime=t_start, location=observer.location)
        
        for name, obj in unique_objects.items():
            alt = obj['coord'].transform_to(frame_start).alt
            if alt.deg > horizon_limit:
                try:
                    target = FixedTarget(coord=obj['coord'], name=name)
                    next_set = observer.target_set_time(t_start, target, which='next', horizon=horizon_limit*u.deg)
                    e_time = min(next_set, t_end) if next_set else t_end
                    
                    if (e_time - t_start).to(u.minute).value > 10:
                        s_dt = t_start.to_datetime(timezone=tz)
                        e_dt = e_time.to_datetime(timezone=tz)
                        
                        results.append({
                            'name': obj['name'],
                            'type': obj['type'],
                            'ra_str': obj['coord'].ra.to_string(unit=u.hour, sep='hms', precision=1),
                            'dec_str': obj['coord'].dec.to_string(unit=u.deg, sep='dms', precision=0, alwayssign=True),
                            'mag': obj['mag'],
                            'observable_str': f"{s_dt.strftime('%H:%M')} - {e_dt.strftime('%H:%M')}"
                        })
                except:
                    pass # Skip if calculation fails
        
        results.sort(key=lambda x: float(x['mag']) if x['mag'] != '-' else 99)
        print(f"Returning {len(results)} observable objects.")
        
        return results, None

    except Exception as e:
        print(f"Sky Query Error: {e}")
        return None, str(e)


# --- Flask Routes ---

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/observable')
def observable():
    return render_template('observable.html')

@app.route('/search_observable', methods=['POST'])
def search_observable():
    lat = float(request.form['latitude'])
    lon = float(request.form['longitude'])
    elev = float(request.form['elevation'])
    start_time_str = request.form['start_time']
    mag_limit = float(request.form['mag_limit'])
    horizon_limit = float(request.form['horizon_limit'])

    results, error = perform_sky_query(lat, lon, elev, start_time_str, mag_limit, horizon_limit)
    
    return render_template('observable.html', results=results, error=error)

@app.route('/deep_dive')
def deep_dive_search():
    return render_template('deep_dive_search.html')

@app.route('/deep_dive_result')
def deep_dive_result():
    object_name = request.args.get('object_name')
    fov_deg = float(request.args.get('fov', 2.0))
    
    if not object_name:
        return "Error: No object name provided.", 400

    object_data, error = get_deep_dive_data(object_name)
    if error:
        return f"Error: {error}", 500

    images = []
    surveys = [
        {'id': 'CDS/P/DSS2/red', 'name': 'DSS2 Red', 'cmap': 'gray'},
        {'id': 'CDS/P/DSS2/blue', 'name': 'DSS2 Blue', 'cmap': 'gray'},
        {'id': 'CDS/P/2MASS/color', 'name': '2MASS Color (IR)', 'cmap': None}, # cmap is ignored for color images
        {'id': 'CDS/P/HI4PI/NH', 'name': 'HI4PI (Hydrogen)', 'cmap': 'inferno'},
        {'id': 'CDS/P/AKARI/FIS/N160', 'name': 'AKARI (Dust)', 'cmap': 'inferno'}
    ]

    for survey in surveys:
        try:
            result = hips2fits.query(
                hips=survey['id'], width=300, height=300,
                ra=object_data['coord'].ra, dec=object_data['coord'].dec,
                fov=fov_deg*u.deg, projection='TAN', format='fits'
            )
            if result:
                img_data = result[0].data
                img_plot = generate_image_plot(img_data, survey['name'], cmap=survey['cmap'])
                images.append({'name': survey['name'], 'plot': img_plot})
        except Exception as e:
            print(f"Could not fetch {survey['name']}: {e}")

    return render_template('deep_dive.html', object=object_data, images=images)

if __name__ == '__main__':
    app.run(debug=True)