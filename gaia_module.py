from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u
import warnings
import numpy as np

# Suppress warnings about VO tables
warnings.filterwarnings("ignore", module="astropy.io.votable.*")

def query_gaia_cone(ra_deg, dec_deg, radius_deg=0.1, limit=50):
    """
    Queries the Gaia DR3 catalog for sources within a cone search.

    Args:
        ra_deg (float): Right Ascension in degrees.
        dec_deg (float): Declination in degrees.
        radius_deg (float): Search radius in degrees (default: 0.1).
        limit (int): Maximum number of rows to return (default: 50).

    Returns:
        list: A list of dictionaries containing source data, or None if query fails.
    """
    coord = SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.degree, u.degree), frame='icrs')
    
    try:
        # Define the query
        # We select standard astrometric parameters and photometry
        query = f"""
        SELECT TOP {limit}
            source_id, ra, dec, parallax, pmra, pmdec, phot_g_mean_mag, bp_rp
        FROM gaiadr3.gaia_source
        WHERE 1=CONTAINS(
            POINT('ICRS', ra, dec),
            CIRCLE('ICRS', {ra_deg}, {dec_deg}, {radius_deg})
        )
        ORDER BY phot_g_mean_mag ASC
        """
        
        # Launch the query asynchronously to avoid blocking if we were in a GUI, 
        # though for Flask sync is usually fine for simple queries. 
        # Using launch_job_async is standard for ADQL queries in astroquery.
        job = Gaia.launch_job_async(query)
        table = job.get_results()
        
        if not table or len(table) == 0:
            return []

        # Convert Astropy Table to list of dicts for easy JSON serialization
        results = []
        for row in table:
            # Handle masked values (replace with None or suitable default)
            def get_val(key):
                if key not in row.colnames:
                    return None
                val = row[key]
                if np.ma.is_masked(val):
                    return None
                return float(val) # Convert numpy types to native python float
            
            results.append({
                'source_id': str(row['source_id']),
                'ra': get_val('ra'),
                'dec': get_val('dec'),
                'parallax': get_val('parallax'),
                'pmra': get_val('pmra'),
                'pmdec': get_val('pmdec'),
                'mag_g': get_val('phot_g_mean_mag'),
                'bp_rp': get_val('bp_rp') # Color index
            })
            
        return results

    except Exception as e:
        print(f"Gaia Query Error: {e}")
        return None

if __name__ == "__main__":
    # Test the function
    print("Testing Gaia Query for Sirius...")
    # Sirius coords: RA 101.287, Dec -16.716
    data = query_gaia_cone(101.287, -16.716, radius_deg=0.05, limit=5)
    if data:
        for star in data:
            print(f"ID: {star['source_id']}, Mag: {star['mag_g']}, Plx: {star['parallax']}")
    else:
        print("No data found or error occurred.")