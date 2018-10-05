"""
By James Simone. Declare the underlying table structure of our databases.
"""
import sys, os
import bz2
import StringIO
import hashlib
import json
import yaml
import datetime

from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.sql import text

Declare = declarative_base()

class DBtimestamp(Declare):
    "DB last modified timestamp"
    __tablename__ = 'modify_times'
    id = Column(Integer, Sequence('modify_time_id_seq'), primary_key=True)
    time = Column(DATETIME)
    def __init__(self,dtobj):
        "UT time now: dtobj = datetime.datetime.utcnow()"
        self.time = dtobj
        return
    pass

class ParameterSet(Declare):
    __tablename__ = 'parameters'
    id = Column(Integer, Sequence('parameter_id_seq'), primary_key=True)
    param = Column(String)
    def __init__(self,param):
        self.param = param
        self.hash()
        return
    def hash(self):
        self.digest = hashlib.md5(self.param).hexdigest()
        return
    def restore(self):
        self.hash()
    def __repr__(self):
        self.restore()
        return '<ParameterSet(%s,"%s")>' % (str(self.id),self.digest)
    pass

def rotateList(v,n):
    "rotate elements of a list"
    #v = v[:] comment out to do in-place
    n = n % len(v)
    head = v[:n]
    v[:n] = []
    v.extend(head)
    return v

class Datum(Declare):
    __tablename__ = 'data'
    id = Column(Integer, Sequence('datum_id_seq'), primary_key=True)
    correlator_id = Column(Integer, ForeignKey('correlators.id'))
    series = Column(String)
    trajectory = Column(Integer)
    tsrc = Column(Integer)
    parameter_id = Column(Integer, ForeignKey('parameters.id'))
    dataBZ2 = Column(BLOB)
    def __init__(self,correlator,series,trajectory,tsrc,params,sdata,doTranslate):
        self.correlator_id = correlator.id
        self.series = series
        self.trajectory = trajectory
        self.tsrc = int(tsrc)
        self.parameter_id = params.id
        if doTranslate:
            rotateList(sdata,tsrc)
            pass
        self.data = '\n'.join(sdata) # newline separated list encoding
        self.dataBZ2 = bz2.compress(self.data,9)
        return
    def __repr__(self):
        sz = len(self.dataBZ2)
        s = '<Datum(' + ','.join(map(str,(self.id,
                                          self.correlator_id,
                                          self.series,
                                          self.trajectory,
                                          self.parameter_id,
                                          len(self.dataBZ2)))) + ')>'
        return (s)
    pass


class Correlator(Declare):
    __tablename__ = 'correlators'
    id = Column(Integer, Sequence('correlator_id_seq'), primary_key=True)
    name = Column(String)
    parameter_id = Column(Integer, ForeignKey('parameters.id'))
    parameters = relationship(ParameterSet,primaryjoin=parameter_id==ParameterSet.id)
    data = relationship(Datum,primaryjoin=id==Datum.correlator_id)
    def __init__(self,name,pset):
        self.name = name
        self.pset = pset
        self.parameter_id = pset.id
        return
    pass

# additional table to support incremental DB builds
class InputFile(Declare):
    __tablename__ = 'correlator_files'
    id = Column(Integer, Sequence('correlator_files_seq'), primary_key=True)
    path = Column(String)
    st_size = Column(Integer)
    st_mtime = Column(Float)
    addtime = Column(DATETIME)
    md5 = Column(String)
    def __init__(self,path,md5=None):
        self.path = path
        st = os.stat(path)
        self.st_size = st.st_size
        self.st_mtime = st.st_mtime # defaults to a float in Python
        self.addtime = datetime.datetime.utcnow()
        self.md5 = md5
        if md5 is None:
            self.md5 = hashlib.md5(open(path,'rb').read()).hexdigest()
        return
    pass
