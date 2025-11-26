#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
MRI DICOM 批量转换为 NIfTI 工具
按原始DICOM文件夹结构组织输出，同一文件夹内的序列放在一起
"""

import os
import sys
import argparse
import logging
import subprocess
import shutil
import tempfile
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Set
from dataclasses import dataclass, field

import pydicom

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class DicomFolder:
    """DICOM文件夹信息"""
    folder_path: str
    folder_name: str
    patient_id: str = ""
    study_date: str = ""
    study_description: str = ""
    series_count: int = 0
    file_count: int = 0
    series_descriptions: List[str] = field(default_factory=list)


class DicomToNiftiConverter:
    """DICOM到NIfTI转换器 - 按文件夹组织输出"""

    def __init__(self,
                 input_dir: str,
                 output_dir: str,
                 compress: bool = True,
                 use_dcm2niix: bool = True,
                 organize_by_folder: bool = True,
                 include_patient_info: bool = True):
        """
        初始化转换器

        Args:
            input_dir: 输入DICOM文件目录
            output_dir: 输出NIfTI文件目录
            compress: 是否压缩输出文件(.nii.gz)
            use_dcm2niix: 是否使用dcm2niix
            organize_by_folder: 按原始文件夹结构组织输出
            include_patient_info: 输出文件名中包含患者信息
        """
        self.input_dir = Path(input_dir).resolve()
        self.output_dir = Path(output_dir).resolve()
        self.compress = compress
        self.use_dcm2niix = use_dcm2niix
        self.organize_by_folder = organize_by_folder
        self.include_patient_info = include_patient_info

        # 存储发现的DICOM文件夹
        self.dicom_folders: Dict[str, DicomFolder] = {}

        if not self.input_dir.exists():
            raise ValueError(f"输入目录不存在: {self.input_dir}")

        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 检查dcm2niix
        if self.use_dcm2niix:
            self._check_dcm2niix()

    def _check_dcm2niix(self) -> bool:
        """检查dcm2niix是否可用"""
        try:
            result = subprocess.run(
                ['dcm2niix', '-v'],
                capture_output=True,
                text=True
            )
            logger.info(f"dcm2niix 可用")
            return True
        except FileNotFoundError:
            logger.warning("dcm2niix 未安装，请先安装：")
            logger.warning("  Ubuntu/Debian: sudo apt-get install dcm2niix")
            logger.warning("  macOS: brew install dcm2niix")
            logger.warning("  Windows: https://github.com/rordenlab/dcm2niix/releases")
            self.use_dcm2niix = False
            return False

    def scan_dicom_folders(self) -> Dict[str, DicomFolder]:
        """
        扫描输入目录，找出所有包含DICOM文件的文件夹
        """
        logger.info(f"扫描DICOM文件夹: {self.input_dir}")

        self.dicom_folders.clear()

        # 遍历所有子目录
        for root, dirs, files in os.walk(self.input_dir):
            dicom_files = []

            # 检查当前目录下的文件
            for file in files:
                file_path = os.path.join(root, file)
                if self._is_dicom_file(file_path):
                    dicom_files.append(file_path)

            # 如果找到DICOM文件，记录这个文件夹
            if dicom_files:
                folder_info = self._analyze_dicom_folder(root, dicom_files)
                if folder_info:
                    self.dicom_folders[root] = folder_info

        # 打印扫描结果
        logger.info(f"\n发现 {len(self.dicom_folders)} 个DICOM文件夹:")
        for folder_path, info in self.dicom_folders.items():
            rel_path = os.path.relpath(folder_path, self.input_dir)
            logger.info(f"  📁 {rel_path}")
            logger.info(f"     患者: {info.patient_id}, 日期: {info.study_date}")
            logger.info(f"     文件数: {info.file_count}, 序列数: {info.series_count}")
            if info.series_descriptions:
                for desc in info.series_descriptions[:5]:  # 最多显示5个
                    logger.info(f"       - {desc}")
                if len(info.series_descriptions) > 5:
                    logger.info(f"       ... 还有 {len(info.series_descriptions) - 5} 个序列")

        return self.dicom_folders

    def _is_dicom_file(self, file_path: str) -> bool:
        """检查文件是否为DICOM文件"""
        try:
            # 快速检查：读取文件头
            with open(file_path, 'rb') as f:
                # 跳过128字节前导
                f.seek(128)
                magic = f.read(4)
                if magic == b'DICM':
                    return True

            # 尝试用pydicom读取
            dcm = pydicom.dcmread(file_path, stop_before_pixels=True, force=True)
            return hasattr(dcm, 'SeriesInstanceUID')

        except Exception:
            return False

    def _analyze_dicom_folder(self, folder_path: str,
                              dicom_files: List[str]) -> Optional[DicomFolder]:
        """分析DICOM文件夹，获取元信息"""
        try:
            # 读取第一个文件获取基本信息
            dcm = pydicom.dcmread(dicom_files[0], stop_before_pixels=True, force=True)

            # 收集所有序列描述
            series_set: Set[str] = set()
            series_uids: Set[str] = set()

            for file_path in dicom_files[:100]:  # 只检查前100个文件以提高速度
                try:
                    d = pydicom.dcmread(file_path, stop_before_pixels=True, force=True)
                    if hasattr(d, 'SeriesInstanceUID'):
                        series_uids.add(str(d.SeriesInstanceUID))
                    if hasattr(d, 'SeriesDescription'):
                        desc = str(d.SeriesDescription).strip()
                        if desc:
                            series_set.add(desc)
                except Exception:
                    continue

            # 创建文件夹信息
            folder_name = os.path.basename(folder_path)

            info = DicomFolder(
                folder_path=folder_path,
                folder_name=folder_name,
                patient_id=str(getattr(dcm, 'PatientID', 'Unknown')).strip(),
                study_date=str(getattr(dcm, 'StudyDate', '')).strip(),
                study_description=str(getattr(dcm, 'StudyDescription', '')).strip(),
                series_count=len(series_uids) if series_uids else len(series_set),
                file_count=len(dicom_files),
                series_descriptions=sorted(list(series_set))
            )

            return info

        except Exception as e:
            logger.warning(f"分析文件夹失败 {folder_path}: {e}")
            return None

    def convert_all(self) -> Dict[str, List[str]]:
        """
        转换所有DICOM文件夹

        Returns:
            字典：{输出文件夹路径: [生成的NIfTI文件列表]}
        """
        if not self.dicom_folders:
            self.scan_dicom_folders()

        if not self.dicom_folders:
            logger.warning("没有找到有效的DICOM文件夹")
            return {}

        results = {}
        total = len(self.dicom_folders)

        for idx, (folder_path, folder_info) in enumerate(self.dicom_folders.items(), 1):
            logger.info(f"\n[{idx}/{total}] 转换: {folder_info.folder_name}")

            try:
                output_subdir, converted_files = self._convert_folder(folder_info)
                results[output_subdir] = converted_files
                logger.info(f"  ✓ 成功生成 {len(converted_files)} 个NIfTI文件")

            except Exception as e:
                logger.error(f"  ✗ 转换失败: {e}")
                results[folder_path] = []

        # 打印总结
        self._print_summary(results)

        return results

    def _convert_folder(self, folder_info: DicomFolder) -> tuple:
        """
        转换单个DICOM文件夹

        Returns:
            (输出目录路径, 生成的文件列表)
        """
        # 确定输出目录结构
        output_subdir = self._get_output_subdir(folder_info)
        output_subdir.mkdir(parents=True, exist_ok=True)

        if self.use_dcm2niix:
            converted_files = self._convert_with_dcm2niix(
                folder_info.folder_path,
                str(output_subdir),
                folder_info
            )
        else:
            converted_files = self._convert_with_python(
                folder_info.folder_path,
                str(output_subdir),
                folder_info
            )

        return str(output_subdir), converted_files

    def _get_output_subdir(self, folder_info: DicomFolder) -> Path:
        """确定输出子目录路径"""
        if self.organize_by_folder:
            # 保持原始目录结构
            rel_path = os.path.relpath(folder_info.folder_path, self.input_dir)

            # 清理路径名
            rel_path = self._clean_path_name(rel_path)

            # 可选：添加患者信息到目录名
            if self.include_patient_info and folder_info.patient_id:
                # 检查路径中是否已包含患者ID
                if folder_info.patient_id not in rel_path:
                    parent = os.path.dirname(rel_path)
                    folder_name = os.path.basename(rel_path)

                    # 构建新的文件夹名：患者ID_日期_原文件夹名
                    parts = []
                    if folder_info.patient_id:
                        parts.append(self._clean_path_name(folder_info.patient_id))
                    if folder_info.study_date:
                        parts.append(folder_info.study_date)
                    parts.append(folder_name)

                    new_folder_name = "_".join(parts)
                    rel_path = os.path.join(parent, new_folder_name) if parent else new_folder_name

            return self.output_dir / rel_path
        else:
            # 扁平结构：所有文件放在输出根目录
            return self.output_dir

    def _clean_path_name(self, name: str) -> str:
        """清理路径/文件名中的非法字符"""
        invalid_chars = '<>:"|?*'
        for char in invalid_chars:
            name = name.replace(char, '_')
        # 移除多余空格
        name = ' '.join(name.split())
        return name.strip()

    def _convert_with_dcm2niix(self, input_folder: str,
                               output_folder: str,
                               folder_info: DicomFolder) -> List[str]:
        """使用dcm2niix转换"""

        # 构建文件名格式
        # %p = 协议名, %s = 序列号, %d = 描述, %n = 患者名
        if self.include_patient_info:
            filename_format = "%p_%s_%d"
        else:
            filename_format = "%s_%d"

        # dcm2niix命令
        cmd = [
            'dcm2niix',
            '-z', 'y' if self.compress else 'n',  # 压缩
            '-f', filename_format,  # 文件名格式
            '-o', output_folder,  # 输出目录
            '-b', 'y',  # 生成bids json
            input_folder  # 输入目录
        ]

        logger.debug(f"执行命令: {' '.join(cmd)}")

        # 执行转换
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            logger.warning(f"dcm2niix警告: {result.stderr}")

        # 收集生成的文件
        output_path = Path(output_folder)
        nifti_files = list(output_path.glob('*.nii.gz' if self.compress else '*.nii'))

        return [str(f) for f in nifti_files]

    def _convert_with_python(self, input_folder: str,
                             output_folder: str,
                             folder_info: DicomFolder) -> List[str]:
        """使用Python库转换（备用方法）"""
        try:
            import SimpleITK as sitk
        except ImportError:
            logger.error("SimpleITK未安装，请运行: pip install SimpleITK")
            return []

        converted_files = []

        # 按序列分组文件
        series_files = defaultdict(list)

        for root, dirs, files in os.walk(input_folder):
            for file in files:
                file_path = os.path.join(root, file)
                try:
                    dcm = pydicom.dcmread(file_path, stop_before_pixels=True, force=True)
                    if hasattr(dcm, 'SeriesInstanceUID'):
                        series_uid = str(dcm.SeriesInstanceUID)
                        series_files[series_uid].append({
                            'path': file_path,
                            'dcm': dcm
                        })
                except Exception:
                    continue

        # 转换每个序列
        for series_uid, files_info in series_files.items():
            try:
                # 获取序列信息
                dcm = files_info[0]['dcm']
                series_desc = str(getattr(dcm, 'SeriesDescription', '')).strip()
                series_num = getattr(dcm, 'SeriesNumber', 0)

                # 生成输出文件名
                name_parts = []
                if series_num:
                    name_parts.append(f"S{series_num:03d}")
                if series_desc:
                    name_parts.append(self._clean_path_name(series_desc))
                else:
                    name_parts.append(series_uid[:8])

                output_name = "_".join(name_parts)
                ext = ".nii.gz" if self.compress else ".nii"
                output_path = os.path.join(output_folder, output_name + ext)

                # 避免重名
                counter = 1
                base_path = output_path.replace(ext, '')
                while os.path.exists(output_path):
                    output_path = f"{base_path}_{counter}{ext}"
                    counter += 1

                # 使用SimpleITK转换
                file_paths = [f['path'] for f in files_info]
                reader = sitk.ImageSeriesReader()
                reader.SetFileNames(sorted(file_paths))
                reader.MetaDataDictionaryArrayUpdateOn()
                reader.LoadPrivateTagsOn()

                image = reader.Execute()
                sitk.WriteImage(image, output_path)

                converted_files.append(output_path)

            except Exception as e:
                logger.warning(f"转换序列 {series_uid[:8]} 失败: {e}")
                continue

        return converted_files

    def _print_summary(self, results: Dict[str, List[str]]):
        """打印转换结果总结"""
        print("\n" + "=" * 70)
        print("转换完成总结")
        print("=" * 70)

        total_folders = len(results)
        total_files = sum(len(files) for files in results.values())
        success_folders = sum(1 for files in results.values() if files)

        print(f"\n📊 统计:")
        print(f"   DICOM文件夹: {total_folders}")
        print(f"   成功转换: {success_folders}")
        print(f"   生成NIfTI文件: {total_files}")

        print(f"\n📁 输出目录: {self.output_dir}")
        print("\n📋 详细结果:")

        for output_folder, files in results.items():
            rel_path = os.path.relpath(output_folder, self.output_dir)
            if rel_path == '.':
                rel_path = '(根目录)'

            if files:
                print(f"\n   ✓ {rel_path}/")
                for f in files:
                    filename = os.path.basename(f)
                    print(f"      - {filename}")
            else:
                print(f"\n   ✗ {rel_path}/ (转换失败)")

        print("\n" + "=" * 70)


def print_folder_structure(converter: DicomToNiftiConverter):
    """打印发现的DICOM文件夹结构"""
    print("\n" + "=" * 70)
    print("DICOM 文件夹结构")
    print("=" * 70)

    for folder_path, info in converter.dicom_folders.items():
        rel_path = os.path.relpath(folder_path, converter.input_dir)
        print(f"\n📁 {rel_path}")
        print(f"   ├── 患者ID: {info.patient_id}")
        print(f"   ├── 检查日期: {info.study_date}")
        print(f"   ├── 检查描述: {info.study_description}")
        print(f"   ├── 文件数量: {info.file_count}")
        print(f"   └── 序列 ({info.series_count}个):")

        for i, desc in enumerate(info.series_descriptions):
            prefix = "       ├──" if i < len(info.series_descriptions) - 1 else "       └──"
            print(f"{prefix} {desc}")

    print("\n" + "=" * 70)


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='批量将MRI DICOM文件转换为NIfTI格式（按文件夹组织输出）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  # 基本使用（推荐）
  %(prog)s -i /path/to/dicom -o /path/to/nifti

  # 查看文件夹结构
  %(prog)s -i /path/to/dicom -o /path/to/nifti --info-only

  # 不包含患者信息
  %(prog)s -i /path/to/dicom -o /path/to/nifti --no-patient-info

  # 扁平输出（所有文件放在同一目录）
  %(prog)s -i /path/to/dicom -o /path/to/nifti --flat

输出目录结构示例:
  output/
  ├── Patient001_20230101_Scan1/
  │   ├── T1_MPRAGE.nii.gz
  │   ├── T2_FLAIR.nii.gz
  │   └── DWI.nii.gz
  └── Patient002_20230102_Scan2/
      ├── T1_MPRAGE.nii.gz
      └── T2_FLAIR.nii.gz
        '''
    )

    parser.add_argument('-i', '--input', required=True,
                        help='输入DICOM文件目录')
    parser.add_argument('-o', '--output', required=True,
                        help='输出NIfTI文件目录')
    parser.add_argument('--no-compress', action='store_true',
                        help='不压缩输出文件（输出.nii而不是.nii.gz）')
    parser.add_argument('--flat', action='store_true',
                        help='扁平输出，不按文件夹组织')
    parser.add_argument('--no-patient-info', action='store_true',
                        help='输出目录名不包含患者信息')
    parser.add_argument('--no-dcm2niix', action='store_true',
                        help='不使用dcm2niix，使用Python库转换')
    parser.add_argument('--info-only', action='store_true',
                        help='仅显示文件夹信息，不进行转换')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='显示详细日志')

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        # 创建转换器
        converter = DicomToNiftiConverter(
            input_dir=args.input,
            output_dir=args.output,
            compress=not args.no_compress,
            use_dcm2niix=not args.no_dcm2niix,
            organize_by_folder=not args.flat,
            include_patient_info=not args.no_patient_info
        )

        # 扫描文件夹
        converter.scan_dicom_folders()

        if not converter.dicom_folders:
            logger.error("没有找到有效的DICOM文件夹")
            sys.exit(1)

        # 仅显示信息
        if args.info_only:
            print_folder_structure(converter)
            sys.exit(0)

        # 显示结构并确认
        print_folder_structure(converter)

        # 执行转换
        print("\n开始转换...\n")
        results = converter.convert_all()

    except KeyboardInterrupt:
        print("\n\n用户取消操作")
        sys.exit(1)
    except Exception as e:
        logger.error(f"程序执行出错: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()