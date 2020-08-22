#! /usr/bin/env python2.7
from __future__ import division, print_function
import sys
import argparse
import dropbox
import os

""" 
Uses the Dropbox python API to download files from shared folders. Avoids needing to 
download everything to your personal computer and then re-upload to a remote server.

Prerequisites:
install dropbox Python package (pip install dropbox)

Create an app: https://www.dropbox.com/developers/apps
Request full access to your user account
Generate an access token on the resulting page

Copy the token and paste it into a file called ~/.dropbox_token
"""

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output_dir", "-o", help="Where to put the downloaded files", required=True)
    inopts = parser.add_mutually_exclusive_group(required=False)
    inopts.add_argument("--url", "-u", help="(OPTIONAL) URL of file to download", required=False, default=None)
    inopts.add_argument("--dir", "-d", help="(OPTIONAL) URL of directory from which to download all files", required=False, default=None)
    parsed = parser.parse_args()
    if parsed.output_dir[-1] == '/':
        parsed.output_dir = parsed.output_dir[0:-1]
    return parsed

def download_dir_contents(dbx, dir, outdir, LIMIT=1000):
    folder_url = dropbox.files.SharedLink(url=dir)
    files = dbx.files_list_folder(path="", shared_link=folder_url, limit=LIMIT)
    base_url = dir.split('?')[0]
    more_files = True
    while more_files:
        for file in files.entries:
            print("Downloading {}...".format(file.name))
            file_url = dropbox.files.SharedLink(url=base_url + "?preview=" + file.name)
            dbx.sharing_get_shared_link_file_to_file(outdir + "/" + file.name, \
                dir, path='/' + file.name)
        if hasattr(files, 'cursor') and files.cursor is not None and len(files.entries) > LIMIT:
            files = dbx.files_list_folder_continue(files.cursor)
        else:
            more_files = False 
            
def main(args):
    """ Main method """
    options = parse_args()
    
    if not os.path.exists(os.getenv('HOME') + '/.dropbox_token'):
        print("ERROR: dropbox token not found.")
        print("Go to https://www.dropbox.com/developers/apps", file=sys.stderr)
        print("If you have not created an app, create one. Request full access to your user account.", file=sys.stderr)
        print("Then go to the app page and generate an access token.", file=sys.stderr)
        print("Create a file in your home directory on the server where you plan to download the file(s) \
called .dropbox_token.", file=sys.stderr)
        print("Copy and paste the access token (and nothing else) into this file.", file=sys.stderr)
        print("For example, echo \"<paste the token here>\" > ~/.dropbox_token", file=sys.stderr)
        exit(1)
        
    token = None
    f = open(os.getenv('HOME') + '/.dropbox_token', 'r')
    for line in f:
        line = line.rstrip()
        if line != '':
            token = line
            break
    f.close()
    
    dbx = dropbox.Dropbox(token)
    
    # If the user has specified a specific file URL, download that. Otherwise, 
    # search through all files/folders shared with the user.
    if options.url is not None:
        file_url = dropbox.files.SharedLink(url=options.url)
        fn = options.url.split('?')[0]
        fn = fn.split('/')[-1]
        print("Downloading {}...".format(fn))
        dbx.sharing_get_shared_link_file_to_file(options.output_dir + '/' + fn, \
            options.url)
    elif options.dir is not None:
        # If the user has provided a link to a dropbox directory instead of a file,
        # attempt to download its contents.
        download_dir_contents(dbx, options.dir, options.output_dir)
    else:
    
        # dbx.sharing_list_folders()
        # dbx.files_list_folder(path)
    
        LIMIT=1000
    
        more_folders = True
        folders = dbx.sharing_list_folders(limit=LIMIT)
    
        while more_folders:
            for folder in folders.entries:
                print("Download contents of {}? (y/n)".format(folder.name))
                user_input = None
                download = False
                while user_input not in ['y', 'Y', 'n', 'N']:
                    user_input = raw_input()
                    if user_input == 'y' or user_input == 'Y':
                        download = True
                    elif user_input == 'n' or user_input == 'N':
                        download = False
                if download:
                    if not os.path.isdir("{}/{}".format(options.output_dir, folder.name)):
                        os.mkdir("{}/{}".format(options.output_dir, folder.name))
                    download_dir_contents(dbx, folder.preview_url, \
                        "{}/{}".format(options, output_dir, folder.name), LIMIT=LIMIT)
        
            if hasattr(folder, 'cursor') and folders.cursor is not None and len(folders.entries) > LIMIT:
                folders = dbx.sharing_list_folders_continue(folders.cursor)
            else:
                more_folders = False

if __name__ == '__main__':
    sys.exit(main(sys.argv))