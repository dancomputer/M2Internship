import requests
import os
import time
import zipfile
from urllib.parse import urlparse, urlencode, urljoin
from netrc import netrc, NetrcParseError
from bs4 import BeautifulSoup

# --- Configuration ---
USERNAME = "daniel54321"  # Replace
PASSWORD = "Hanakuso007!"  # Replace
OUTPUT_DIR = "gpw_africa_data" 
SEDAC_BASE_URL = "http://sedac.ciesin.columbia.edu/gpw/app/"
DELAY_SECONDS = 1

# --- 1. Earthdata Login (.netrc setup) ---

def setup_netrc(username, password):
    """Creates or updates the .netrc file with Earthdata Login credentials."""
    netrc_path = os.path.expanduser("~/.netrc")
    machine = 'urs.earthdata.nasa.gov'

    try:
        with open(netrc_path, 'r') as f:
            content = f.read()
        if f'machine {machine}' in content:
            print(f".netrc file already contains Earthdata Login credentials at: {netrc_path}")
            return
        with open(netrc_path, 'a') as f:
            f.write(f'\nmachine {machine}\n\tlogin {username}\n\tpassword {password}\n')
        os.chmod(netrc_path, 0o600)
        print(f".netrc file updated (appended) at: {netrc_path}")
    except FileNotFoundError:
        with open(netrc_path, 'w') as f:
            f.write(f'machine {machine}\n\tlogin {username}\n\tpassword {password}\n')
        os.chmod(netrc_path, 0o600)
        print(f".netrc file created at: {netrc_path}")

# --- 2. SEDAC Download URL Generation (Corrected for Form Submission) ---
def get_sedac_download_url(session, base_url, params):
    """
    Constructs the SEDAC application URL, makes an initial GET request,
    then simulates a POST request (form submission) to get the download URL.
    """
    # 1. Initial GET request (with parameters).
    initial_url = base_url + "?" + urlencode(params)
    print(f"Initial GET request: {initial_url}")
    initial_response = session.get(initial_url, allow_redirects=True)
    initial_response.raise_for_status()

    # 2. Prepare POST data (mimic form submission).  Crucially,
    #    we get the form data from the *response* to the GET request.
    soup = BeautifulSoup(initial_response.text, 'html.parser')
    form = soup.find('form')  # Find the form
    if not form:
      raise Exception("Form not found")

    post_data = {}
    for input_tag in form.find_all('input'):
        if input_tag.get('name'): # Ensure the input has a name.
            post_data[input_tag['name']] = input_tag.get('value', '')
    for select_tag in form.find_all('select'):  # Include select tags.
        selected_option = select_tag.find('option', selected=True)
        post_data[select_tag['name']] = selected_option.get('value', '') if selected_option else ''

     # Add the 'data' (grid type), resolution, and year from our *params*.
    post_data['data'] = params['data']
    post_data['resolut'] = params['resolut']
    post_data['year'] = params['year']
    #SEDAC also requires this param
    post_data['submit'] = 'submit'


    # 3. Make the POST request (to the same URL as the form's action).
    post_url = urljoin(base_url, form['action'])  # Get form's action URL
    print(f"POST request: {post_url} with data: {post_data}")
    post_response = session.post(post_url, data=post_data, allow_redirects=True)
    post_response.raise_for_status()

    # 4.  *Now* find the download link in the POST response.
    post_soup = BeautifulSoup(post_response.text, 'html.parser')
    download_link = post_soup.find('a', href=lambda href: href and 'downloads/data' in href)

    if download_link:
        final_url = urljoin(base_url, download_link.get('href'))
        print(f"Final download URL: {final_url}")
        return final_url
    else:
        raise Exception(f"Download link not found after POST request to: {post_response.url}")


# --- 3. Download and Extract ---
def download_and_extract(session, urls, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for url in urls:
        filename = os.path.basename(urlparse(url).path)
        filepath = os.path.join(output_dir, filename)

        if not os.path.exists(filepath):
            print(f"Downloading: {url}")
            try:
                response = session.get(url, stream=True)
                response.raise_for_status()

                with open(filepath, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                print(f"Downloaded: {filename}")

                if filename.endswith('.zip'):
                    print(f"Extracting: {filename}")
                    try:
                        with zipfile.ZipFile(filepath, 'r') as zip_ref:
                            extract_path = os.path.abspath(output_dir)
                            zip_ref.extractall(path=extract_path)
                            print(f"Extracted to: {extract_path}")
                        os.remove(filepath)
                    except zipfile.BadZipFile:
                        print(f"Error: {filename} is not a valid ZIP file.")
                    except Exception as e:
                        print(f"Error extracting {filename}: {e}")
            except requests.exceptions.RequestException as e:
                print(f"Error downloading {url}: {e}")

# --- Main Execution ---
if __name__ == "__main__":
    setup_netrc(USERNAME, PASSWORD)
    session = requests.Session()

    urls = []
    base_params = {
        'region': 'Africa',
        'type': 'tif',
        'level-scope': 'continent',
    }
    
    files = ['gpwv3', 'grumpv1']
    data_types = {  # Map file names to their corresponding 'data' values
        'gpwv3': ['pcount'],#, 'pdens', 'uadle', 'uald', 'landar', 'ctryar'],
        'grumpv1': ['pcount']#, 'pdens', 'uadle']
    }

    resolutions = ['2pt5', '15', '30', '1']
    years = ['1990','1995','2000'] #GRUMP only goes to 2000
    years_gpw = ['1990','1995', '2000', '2005', '2010', '2015', '2020']

    for file in files:
        for data_type in data_types[file]:
            for resolution in resolutions:
                for year in (years_gpw if file == 'gpwv3' else years): # Year loop
                    params = base_params.copy()
                    params['file'] = file
                    params['data'] = data_type
                    params['resolut'] = resolution
                    params['year'] = year
                    try:
                        download_url = get_sedac_download_url(session, SEDAC_BASE_URL, params)
                        urls.append(download_url)
                        time.sleep(DELAY_SECONDS)
                    except Exception as e:
                        print(f"Skipping URL due to error: {e}")
                        continue


    download_and_extract(session, urls, OUTPUT_DIR)
    print("Download and Extraction Complete.")
