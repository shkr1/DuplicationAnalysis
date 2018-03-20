import os
import subprocess
from .grapher import run
from contextlib import contextmanager
from flask import Flask, request, session, g, redirect, url_for, abort, \
        render_template, flash, send_from_directory
from werkzeug.utils import secure_filename


app = Flask(__name__) # create the application instance :)
app.config.from_object(__name__) # load config from this file

UPLOAD_FOLDER = os.path.join(app.root_path, 'duplications')
RESULT_FOLDER = os.path.join(app.root_path, 'duplications/Result')
ALLOWED_EXTENSION = ['fasta']

# Load default config and override config from an environment variable
app.config.update(dict(
    SECRET_KEY='development key',
    USERNAME='admin',
    PASSWORD='default',
    UPLOAD_FOLDER=UPLOAD_FOLDER,
    RESULT_FOLDER=RESULT_FOLDER
))
app.config.from_envvar('DUPLICATION_SETTINGS', silent=True)

@app.route('/')
def index():
    return render_template('index.html')

def allowed_file(filename):
    return '.' in filename and \
            filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSION

def apply_duplication(filename, gen):
    prevdir = os.getcwd()
    os.chdir(UPLOAD_FOLDER)
    try:
        script = "run_duplication_analysis -i" + filename
        with open("GET_GENES.txt", "w") as f:
            f.write(gen)
        subprocess.call("bash " + script, shell=True)
        run(RESULT_FOLDER + '/' + "dS_" + filename + "_RESULT.txt")

    finally:
        os.chdir(prevdir)

@app.route('/add', methods=['POST'])
def upload_entry():
    # Check if the post request has the file part
    if len(request.files) > 1:
        flash('Error inesperado')
        return redirect(url_for('index'))

    if 'file' not in request.files:
        flash('No file part')
        return redirect(url_for('index'))

    # file = request.files['file']
    
    # if user does not select file, browser also
    # submit an empty part without filename
    for file in request.files.getlist('file'):
        if file.filename == '':
            flash('No selected file')
            return redirect(url_for('index'))

        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        else:
            ext = file.filename.split('.')[-1]
            flash('.' + ext + ' extensión no soportada')
            return redirect(url_for('index'))
    name = request.form['name']
    gen = request.form['gen']
    apply_duplication(name, gen)
    flash('uploads/dS_' + name + '_RESULT.png')
    return redirect(url_for('index'))

@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(app.config['RESULT_FOLDER'],
                               filename)