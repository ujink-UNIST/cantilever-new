# cantilever-new

`modal3/uniform3`, `modal4/uniform4`, `modal8/uniform8` 형상을 pygalmesh volume mesh → ANSYS SOLID187 TET10으로 변환하고, 네 재료의 무구속 free-free modal 해석을 실행합니다.

## 입력 형상

OBJ/STL 모델은 용량 때문에 Git 저장소가 아닌 `v1.0` Release에 압축해 두었습니다. `run.py`는 입력 형상이 없으면 다음 파일을 자동 다운로드·압축 해제합니다.

- `modal_cantilever3.obj.zip`
- `uniform_cantilever3.obj.zip`
- `modal_cantilever4.obj.zip`
- `uniform_cantilever4.obj.zip`
- `modal-8.stl.zip`
- `uniform-8.stl.zip`

공개 저장소이므로 GitHub 인증은 필요하지 않습니다. 기존 3번 모델의 최종 OBJ 실체 부피는 modal `83,319.77 mm³`, uniform `83,341.32 mm³`로 차이는 `0.026%`입니다. uniform 공통 반지름은 `0.806104 mm`입니다.

## 8번 STL 모델

8번 모델은 OBJ 변환 없이 Binary STL을 블록 단위로 읽습니다. `models8.json`에 원본/ZIP SHA-256과 형상 정보를 기록했습니다. 로컬 STL은 저장소 루트에 두며 Git에는 포함하지 않습니다.

| 모델 | 삼각형 수 | STL 부피 (mm³) |
|---|---:|---:|
| `modal-8.stl` | 49,735,004 | 85,595.7088 |
| `uniform-8.stl` | 48,574,016 | 85,419.9356 |

두 모델은 지그를 제외한 1:1 mm 형상이며, 입력 부피 차이는 modal 기준 약 0.205%입니다. 프린팅용 지그 포함 부피와 구별해야 합니다. 두 모델 모두 voxel 간격은 0.2 mm, 기존 3/4 모델 설정은 유지합니다. 실제 STL→voxel 검증에서 부피 변화는 modal −0.160%, uniform −0.200%이며 odd ray는 모두 0개였습니다. 이는 tetra 생성/ANSYS 실행 또는 주파수 수렴 검증 결과는 아닙니다.

8번 두 모델만 실행:

```powershell
python .\run.py --cases modal-8 uniform-8
```

실제 해석 없이 8번 whitened 입력 파일만 생성:

```powershell
python .\run_all_modal.py --cases modal-8 uniform-8 --prepare-only
```

전체 모델은 인자 없이 실행합니다. 체적 메시의 기존 허용 오차(voxel 부피 15%, tetra 부피 10%)는 실패 감지 기준이지 정확도 보증이 아닙니다. 주파수 비교에는 두 모델의 실제 체적 메시 부피와 요소 크기 수렴을 별도로 확인해야 합니다.

## 실행

Windows에서는 CGAL 의존성까지 설치되는 conda-forge 환경을 사용합니다.

```powershell
conda env create -f environment.yml
conda activate cantilever
python .\run.py
```

ANSYS를 자동으로 찾지 못하면:

```bash
python run.py --exe "C:\\Program Files\\ANSYS Inc\\v252\\ansys\\bin\\winx64\\ANSYS252.exe"
```

`run.py`가 다음 작업을 순서대로 수행합니다.

1. OBJ/Binary STL을 pygalmesh TET4 volume mesh로 변환
2. edge 중간 절점을 추가해 SOLID187 TET10 CDB 생성
3. 각 형상에서 whitened `ν=0.30`, `ν=0.39` 두 번만 해석
4. `sqrt(E/ρ)` scaling으로 네 재료의 주파수 생성
5. `runs/free_free_modal_frequencies_all.csv` 저장

| 재료 | E (MPa) | ν | ρ (tonne/mm³) |
|---|---:|---:|---:|
| Inconel 718 | 165,000 | 0.30 | 8.220e-9 |
| Structural Steel | 200,000 | 0.30 | 7.850e-9 |
| EOS MaragingSteel MS1 (열처리) | 190,000 | 0.30 | 8.100e-9 |
| Formlabs Nylon 12 | 1,850 | 0.39 | 9.900e-10 |

ANSYS 기준 물성은 `E=1 MPa`, `ρ=1 tonne/mm³`이며 포아송비별로 해석합니다. 같은 포아송비의 실제 주파수는 `f = f_white × sqrt(E/ρ)`로 계산합니다. 경계조건은 없는 true free-free입니다. 추출된 12개 모드 중 1–6차는 강체 모드, 7–12차는 탄성 고유진동 `f1–f6`입니다.

## s0.5 nTop thick graph

두 CSV 그래프를 nTop에서 가져올 수 있는 3MF beam-lattice/ball thick graph로 변환합니다.

```powershell
python .\csv_to_thick_graph_3mf.py
```

출력은 `geometry/*_s0.5_thick_lattice_graph.3mf`입니다. 설계 모델 반지름 범위는 `0.25000–0.55201 mm`, 최종 uniform 공통 반지름은 `0.411625 mm`입니다. 최종 OBJ 실체 부피는 modal `10,658.24 mm³`, uniform `10,657.99 mm³`로 차이는 `0.0023%`입니다.

## 3D 프린팅용 경량 OBJ

Release의 `*_print.obj.zip`은 원본 대비 face 수를 약 90% 줄이고 MeshFix로 폐쇄·manifold 상태를 복구한 파일입니다.

| 형상 | 원본 face | 경량 face | 부피 변화 |
|---|---:|---:|---:|
| modal | 49,907,540 | 4,987,952 | -0.177% |
| uniform | 48,789,056 | 4,875,214 | -0.184% |

두 경량 모델 사이의 부피 차이는 약 `0.019%`입니다. 다시 만들려면 `python .\reduce_obj_for_print.py`를 실행합니다.
