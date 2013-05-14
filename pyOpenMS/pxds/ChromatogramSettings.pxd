from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string
from InstrumentSettings cimport *
from Precursor cimport *
from Peak1D cimport *
from SourceFile cimport *
from PeptideIdentification cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DataValue cimport *
from String cimport *
from Types cimport *
from DataProcessing cimport *

from Product cimport *
from AcquisitionInfo cimport *

cdef extern from "<OpenMS/METADATA/ChromatogramSettings.h>" namespace "OpenMS":

    cdef cppclass ChromatogramSettings:

        ChromatogramSettings()    nogil except +
        ChromatogramSettings(ChromatogramSettings)    nogil except +

        Precursor getPrecursor()
        void setPrecursor(Precursor p)
        Product getProduct()
        void setProduct(Product p)
        String getNativeID()
        void setNativeID(String native_id)

        #libcpp_vector[DataProcessing] getDataProcessing
        #void setDataProcessing(libcpp_vector[DataProcessing data_processing] )
        libcpp_vector[DataProcessing] getDataProcessing() nogil except +
        void setDataProcessing(libcpp_vector[DataProcessing])   nogil except +

