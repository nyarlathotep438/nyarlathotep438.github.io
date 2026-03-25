# 2026年3月25日 研究日志
今天也倒腾了一天数据，给大家分享一下我偶尔会用的办公自动化Python代码吧

当你需要批量生成一些内容相近的参数文件时，你可以直接用这个↓


```python
for i in range(1, 31):
    # 生成文件名
    filename = f"music_{i:02d}.sdef"
    
    # 生成内容
    content = f"""--
wave = "/Effects/Cockpit/DPlayer/music_{i:02d}"
inner_radius = 10
outer_radius = 100"""
    
    # 写入文件
    with open(filename, "w", encoding="utf-8") as f:
        f.write(content)
```

当你需要把某个文件夹的文件批量重命名，你可以用这个↓


```python
import os
import glob

def rename_wav_files(directory):
    # 获取目录中所有wav文件
    wav_files = glob.glob(os.path.join(directory, "*.wav"))
    
    # 按文件修改时间排序（最早的文件为01）
    wav_files.sort(key=os.path.getmtime)
    
    # 如果没有wav文件则退出
    if not wav_files:
        print("目录中没有找到wav文件")
        return
    
    # 重命名所有文件
    for i, file_path in enumerate(wav_files, 1):
        # 生成新文件名
        new_name = f"music_{i:02d}.wav"
        new_path = os.path.join(directory, new_name)
        
        # 重命名文件
        try:
            os.rename(file_path, new_path)
            print(f"重命名成功: {os.path.basename(file_path)} -> {new_name}")
        except Exception as e:
            print(f"重命名失败: {file_path} | 错误: {str(e)}")

if __name__ == "__main__":
    target_dir = r"C:\Users\Admin\Music\Ready2Change"
    
    # 确认目录存在
    if not os.path.exists(target_dir):
        print(f"目录不存在: {target_dir}")
    elif not os.path.isdir(target_dir):
        print(f"路径不是目录: {target_dir}")
    else:
        rename_wav_files(target_dir)
        print("\n所有文件重命名完成！")
```

*今天干的活实在不方便展示，今天就先这样吧，下个博客再见！*
