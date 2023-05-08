import os
import tempfile
import shutil


### build a FIFO to handle temp files
class FIFO:
    def __init__(self, dir, name, prefix="temp_file_dir", suffix = ""):
        self.tempdir = tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=dir)
        self.filename = os.path.join(self.tempdir, name)
    def __enter__(self):
        if os.path.exists(self.filename):
            os.unlink(self.filename)
        os.mknod(self.filename)
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        if hasattr(self, "tempdir"):
            shutil.rmtree(self.tempdir)
