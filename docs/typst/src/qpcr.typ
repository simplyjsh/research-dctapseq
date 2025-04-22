#import "@local/huangjac:1.0.0":*

= qPCR Analysis Tools used in DC-TAP

== Installation

Before we download and run the qPCR analysis script, check if your
developer environment is properly setup.

*Developer Environment Setup*

- We recommend working with the Command Prompt using the #link("https://apps.microsoft.com/detail/9n0dx20hk701?hl=en-US&gl=US")[Terminal App] which you can download off the Microsoft Store.

- We also recommend the #link("https://code.visualstudio.com/download")[Visual Studio Code] as your text editor.

- Ensure that *Python* is installed on your machine. You can check if it installed with the command below:

#codly(header: [Check Python Installation])
```powershell
python -V
```

*Downloading and Installing the Script*

1. Download the qPCR Analysis Script.

Note: This script is still a work in-progress inside a private git repo 
hence the roundabout method of downloading a compressed file. Also, the script
only supports ΔΔCt

#codly(header: [Downloading and Unzipping Script])
```powershell
$zipUrl = "https://drive.google.com/uc?id=1efuTbsrjgHvpBp1zUw_U-RhnmE96hX8_"
$zipPath = "$env:TEMP\tempfile.zip"
$destPath = "$PWD"
Invoke-WebRequest -Uri $zipUrl -OutFile $zipPath; Expand-Archive -Path $zipPath -DestinationPath $destPath -Force; Remove-Item $zipPath
```

2. Setup a Python Virtual Environment

Change your directory to the qPCR script's root directory. Then create a 
virtual environment with the following commands below.

#codly(header: [Create a Python Virtual Environment])
```powershell
python -m venv .venv
. .\.venv\Scripts\activate
```

Once you create the virtual environment, install the required Python
packages listed in _requirements.txt_.

#codly(header: [Setup Virtual Environment])
```powershell
pip install --upgrade pip
pip install -r requirements.txt
```

*Complete*

You can start editing and running the qPCR analysis scripts.

#codly(header: [Open IDE])
```powershell
code .
```
== Running the qPCR Script

*Before you begin*

Every time you exit out of your editor or the CMD or restart the qPCR analysis,
you *MUST* activate the python virtual environment first before running 
the script.

*1. Structuring the data*

In order to process your data, we need to know where the data lives and be
organized in a format that the script expects the data to be in.

The complete qPCR data files we required are listed below:

1. Raw Ct values
2. Template layout
3. Primer layout

The expected structure we expect the data to be organized in are detailed
below:

#set table(
  stroke: none,
  gutter: 0.2em,
  fill: (x, y) =>
    if x == 0 or y == 0 { gray },
  inset: (right: 1.5em),
)

#show table.cell: it => {
  if it.x == 0 or it.y == 0 {
    set text(white)
    strong(it)
  } else if it.body == [] {
    // Replace empty cells with 'N/A'
    pad(..it.inset)[_N/A_]
  } else {
    it
  }
}

1. Raw Ct values Table
#table(
  columns: 4,
  [], [Exam 1], [Exam 2], [Exam 3],

  [John], [], [], [],
  [Mary], [], [], [],
  [Robert], [], [], [],
)
