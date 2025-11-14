"""
    制作CoLM运行的站点数据
        CoLM：站点数据->适用于原本的CROP模型
        CSM：站点数据->适用于新的CROP模型
"""
import os
import numpy as np
import xarray as xr

def Srfdata_CoLMCROP(srfdata_filedir, site, origfile_filepath):
    """
        生成单点的驱动数据
    """
    if len(origfile_filepath) > 0:
        #   1 读取原始文件
        srfdata = xr.open_dataset(origfile_filepath)
        #   2 修改原始文件中的部分信息
        srfdata["latitude"].values = np.array(site["lat"], dtype=np.float32)
        srfdata["longitude"].values = np.array(site["lon"], dtype=np.float32)
        srfdata["pfttyp"].values = np.array([site["crop"]], dtype=np.int32)
        srfdata["pctpfts"].values = np.array([1.0], dtype=np.float32)
        srfdata["croptyp"] = xr.DataArray(np.array([site["crop"] - 14], dtype=np.int32), dims=["pft"])
        srfdata["pctcrop"] = xr.DataArray(np.array([1.0], dtype=np.float32), dims=["pft"])
        #   3 保存数据
        srfdata_filepath = os.path.join(srfdata_filedir, "srfdata_{}_CoLMCROP.nc".format(site["siteid"]))
        srfdata.to_netcdf(srfdata_filepath)
    else:
        #   1 生成输出文件
        srfdata_filepath = os.path.join(srfdata_filedir, "srfdata_{}_CoLMCROP.nc".format(site["siteid"]))
        #   2 数据生成xarray
        srfdata = xr.Dataset({"latitude" : np.array(site["lat"], dtype=np.float32),
                            "longitude": np.array(site["lon"], dtype=np.float32),
                            "USGS_classification":3,
                            "IGBP_classification":12,
                            "pfttyp" : (["pft"], np.array([site["crop"]], dtype=np.int32)),
                            "pctpfts": (["pft"], np.array([1.0], dtype=np.float32)),
                            "croptyp": (["pft"], np.array([site["crop"] - 14], dtype=np.int32)),
                            "pctcrop": (["pft"], np.array([1.0], dtype=np.float32))})
        #   3 保存数据
        srfdata.to_netcdf(srfdata_filepath)

def Srfdata_CSM(srfdata_filedir, site, origfile_filepath):
    """
        生成单点的驱动数据
    """
    if len(origfile_filepath) > 0:
        #   1 读取原始文件
        srfdata = xr.open_dataset(origfile_filepath)
        #   2 修改原始文件中的部分信息
        srfdata["latitude"].values = np.array(site["lat"], dtype=np.float32)
        srfdata["longitude"].values = np.array(site["lon"], dtype=np.float32)
        srfdata["pfttyp"].values = np.array([15], dtype=np.int32)
        srfdata["pctpfts"].values = np.array([1.0], dtype=np.float32)
        srfdata["croptyp"] = xr.DataArray(np.array([site["crop"]], dtype=np.int32), dims=["pft"])
        srfdata["pctcrop"] = xr.DataArray(np.array([1.0], dtype=np.float32), dims=["pft"])
        #   3 保存数据
        srfdata_filepath = os.path.join(srfdata_filedir, "srfdata_{}_CSM.nc".format(site["siteid"]))
        srfdata.to_netcdf(srfdata_filepath)
    else:
        #   1 生成输出文件
        srfdata_filepath = os.path.join(srfdata_filedir, "srfdata_{}_CSM.nc".format(site["siteid"]))
        #   2 数据生成xarray
        srfdata = xr.Dataset({"latitude" : np.array(site["lat"], dtype=np.float32),
                            "longitude": np.array(site["lon"], dtype=np.float32),
                            "USGS_classification":3,
                            "IGBP_classification":12,
                            "pfttyp" : (["pft"], np.array([15], dtype=np.int32)),
                            "pctpfts": (["pft"], np.array([1.0], dtype=np.float32)),
                            "croptyp": (["pft"], np.array([site["crop"]], dtype=np.int32)),
                            "pctcrop": (["pft"], np.array([1.0], dtype=np.float32))})
        #   3 保存数据
        srfdata.to_netcdf(srfdata_filepath)


def Srfdata_CoLM(srfdata_filedir, site, origfile_filepath):
    """
        生成单点的驱动数据
    """
    if len(origfile_filepath) > 0:
        #   1 读取原始文件
        srfdata = xr.open_dataset(origfile_filepath)
        #   2 修改原始文件中的部分信息
        srfdata["latitude"].values = np.array(site["lat"], dtype=np.float32)
        srfdata["longitude"].values = np.array(site["lon"], dtype=np.float32)
        srfdata["pfttyp"].values = np.array([15], dtype=np.int32)
        srfdata["pctpfts"].values = np.array([1.0], dtype=np.float32)
        #   3 保存数据
        srfdata_filepath = os.path.join(srfdata_filedir, "srfdata_{}_CoLM.nc".format(site["siteid"]))
        srfdata.to_netcdf(srfdata_filepath)
    else:
        #   1 生成输出文件
        srfdata_filepath = os.path.join(srfdata_filedir, "srfdata_{}_CoLM.nc".format(site["siteid"]))
        #   2 数据生成xarray
        srfdata = xr.Dataset({"latitude" : np.array(site["lat"], dtype=np.float32),
                            "longitude": np.array(site["lon"], dtype=np.float32),
                            "USGS_classification":3,
                            "IGBP_classification":12,
                            "pfttyp" : (["pft"], np.array([15], dtype=np.int32)),
                            "pctpfts": (["pft"], np.array([1.0], dtype=np.float32))})
        #   3 保存数据
        srfdata.to_netcdf(srfdata_filepath)


if __name__ == "__main__":
    #   1 基础信息
    srfdata_filedir = "/stu01/lianghb21/colm_crop/model/CoLM_CSM/data/"
    if not os.path.exists(srfdata_filedir):
        os.makedirs(srfdata_filedir)
    #   2 站点信息
    lat = 41.1481
    lon = 121.2017
    siteid = "JinZhou"
    croptype = "maize"
    origfile_filepath = ""
    #   2 生成站点数据，给定站点信息 CoLM
    croptype_colm = {"maize": 17, "soybean": 23, "spwheat": 19, "wtwheat": 21, "rice": 61}
    site = {"siteid": siteid, "lat": lat, "lon": lon, "crop": croptype_colm[croptype]}
    Srfdata_CoLM(srfdata_filedir, site, origfile_filepath)
    Srfdata_CoLMCROP(srfdata_filedir, site, origfile_filepath)
    #   3 生成站点数据，给定站点信息 CSM
    croptype_csm = {"maize": 1, "soybean": 3, "spwheat": 2, "wtwheat": 2, "rice": 4}
    site = {"siteid": siteid, "lat": lat, "lon": lon, "crop": croptype_csm[croptype]}
    Srfdata_CSM(srfdata_filedir, site, origfile_filepath)
    