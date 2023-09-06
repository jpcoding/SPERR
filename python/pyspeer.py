import ctypes
import os
import sys  
import numpy as np 
class pysperr:
    def __init__(self, libpath):
        self.sperr = ctypes.cdll.LoadLibrary(libpath)
        # self.sperr.sperr_comp_2d = self.sperr.sperr_comp_2d
        self.sperr.sperr_comp_2d.argtypes = (ctypes.c_void_p, ctypes.c_int32,
                                             ctypes.c_size_t, ctypes.c_size_t, ctypes.c_int32, ctypes.c_double,
                                             ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_size_t))

        self.sperr.sperr_comp_2d.restype = ctypes.c_int32
        
        # self.sperr.sperr_comp_3d = self.sperr.sperr_comp_3d
        self.sperr.sperr_comp_3d.argtypes = [ctypes.c_void_p, ctypes.c_int32,
                            ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t, ctypes.c_int32, ctypes.c_double,
                            ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_size_t)]
        self.sperr.sperr_comp_3d.restype = ctypes.c_int32
        
        # self.sperr.sperr_decomp_2d = self.sperr.sperr_decomp_2d
        self.sperr.sperr_decomp_2d.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int32,
                            ctypes.POINTER(ctypes.c_size_t), ctypes.POINTER(ctypes.c_size_t),
                            ctypes.POINTER(ctypes.c_void_p)]
        self.sperr.sperr_decomp_2d.restype = ctypes.c_int32
        
        # self.sperr.sperr_decomp_3d = self.sperr.sperr_decomp_3d
        self.sperr.sperr_decomp_3d.argtypes =[ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int32,
                            ctypes.POINTER(ctypes.c_size_t), ctypes.POINTER(ctypes.c_size_t),ctypes.POINTER(ctypes.c_size_t),
                            ctypes.POINTER(ctypes.c_void_p)]
        self.sperr.sperr_decomp_3d.restype = ctypes.c_int32
        

    def compress(self, src, mode_str, quality):
        # src has to be a numpy array check at the beginning
        if not isinstance(src, np.ndarray):
            raise ValueError("Input has to be a numpy array")
        dst = ctypes.c_void_p()
        dst_len = ctypes.c_size_t()
        modes = {"bpp": 1, "psnr": 2, "pwe": 3}
        if mode_str not in modes:
            raise ValueError("Invalid mode")
        mode = modes[mode_str]
        is_float = 1 if src.dtype == np.float32 else 0
        shape = src.shape
        if len(src.shape) == 2:
            result = self.sperr.sperr_comp_2d(src.ctypes.data, is_float, shape[1], 
                           shape[0], mode, quality, ctypes.byref(dst), ctypes.byref(dst_len))
        elif len(src.shape) == 3:
            result = self.sperr.sperr_comp_3d(src.ctypes.data, is_float, shape[2],shape[1],shape[0], 
                         mode, quality, ctypes.byref(dst), ctypes.byref(dst_len))
        else:
            raise ValueError("Only 2D or 3D arrays are supported")
        if result != 0:
            raise ValueError("Compression failed")
        outbytes = ctypes.create_string_buffer(dst_len.value)
        ctypes.memmove(outbytes, dst, dst_len.value)
        ## free the memory of dst
        self.sperr .free(dst)
        dst = None
        dst_bytes = None
        # outbytes = np.frombuffer(outbytes, dtype=np.uint8)
        print(type(outbytes))
        return outbytes, dst_len.value
    
    def decompress(self,src, is_float, dims):
        src_ptr = ctypes.cast(src, ctypes.c_void_p)
        dtype = np.float32 if is_float == 1 else np.float64
        dst = ctypes.c_void_p()
        dtype = np.float32 if is_float == 1 else np.float64
        size = dims[0] * dims[1] * (dims[2] if len(dims) == 3 else 1)
        if len(dims) == 2:
            result = self.sperr.sperr_decomp_2d(src_ptr, ctypes.c_size_t(len(src)), 
                             ctypes.c_int32(is_float), 
                             ctypes.byref(ctypes.c_size_t(dims[1])),
                             ctypes.byref(ctypes.c_size_t(dims[0])), ctypes.byref(dst))
        elif len(dims) == 3:
            result = self.sperr.sperr_decomp_3d(src_ptr, ctypes.c_size_t(len(src)), 
                             ctypes.c_int32(is_float),  
                             ctypes.byref(ctypes.c_size_t(dims[2])),
                             ctypes.byref(ctypes.c_size_t(dims[1])),
                             ctypes.byref(ctypes.c_size_t(dims[0])), ctypes.byref(dst))
        else:
            raise ValueError("Only 2D or 3D arrays are supported")
        if result != 0:
            raise ValueError("Decompression failed")
        
        outbytes = ctypes.create_string_buffer(size *(4 if is_float else 8))
        ctypes.memmove(outbytes, dst, size *(4 if is_float else 8))
        dst_array = np.frombuffer(outbytes, dtype=dtype, count=size).reshape(dims)
        self.sperr .free(dst)
        src_ptr = None
        dst = None
        outbytes = None
        return dst_array
    
    def compress_decompress(self, src, mode_str, quality):
        dst, dst_len = self.compress(src, mode_str, quality)
        is_float = 1 if src.dtype == np.float32 else 0
        dims = src.shape
        dst_array = self.decompress(dst, is_float, dims)
        result = {}
        result["ddata"] = dst_array
        result["cr"] = 1.0* src.nbytes/dst_len
        return result 
            
            