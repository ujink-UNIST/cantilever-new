# cantilever-new

`modal3/uniform3`, `modal4/uniform4`, `modal8/uniform8` 형상을 pygalmesh volume mesh → ANSYS SOLID187 TET10으로 변환하고, 네 재료의 무구속 free-free modal 해석을 실행합니다.

## 입력 형상

3/4번 모델은 Git에 포함된 3MF beam-lattice/ball graph를 직접 사용합니다. nTop 및 중간 OBJ/surface mesh는 해석에 필요하지 않습니다.

- `geometry/modal_cantilever3_corrected.3mf`
- `geometry/uniform_cantilever3_volume_matched.3mf`
- `geometry/*_s0.5_thick_lattice_graph.3mf`

8번 모델만 용량 때문에 `v1.0` Release의 `modal-8.stl.zip`, `uniform-8.stl.zip`을 자동 다운로드합니다. 공개 저장소이므로 GitHub 인증은 필요하지 않습니다.

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

전체 모델은 인자 없이 실행합니다. 3MF voxel 쌍은 부피 차이 1%, 최종 tetra는 목표 부피 차이 3%를 넘으면 실패합니다. 주파수 비교에는 요소 크기 수렴을 별도로 확인해야 합니다.

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

1. 3MF graph 또는 Binary STL을 voxel union → pygalmesh/CGAL TET4로 변환
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

## 3MF 직접 voxel 부피 검증

3MF의 butt-cap beam과 node ball을 직접 raster union합니다. 중간 surface가 없어 watertight/self-intersection 문제에 의존하지 않습니다.

| 쌍 | Modal voxel (mm³) | Uniform voxel (mm³) | 차이 |
|---|---:|---:|---:|
| 3 | 86,399.552 | 86,428.912 | 0.034% |
| 4 | 10,805.797 | 10,813.515 | 0.071% |

최소 스트럿 지름은 modal에서 4 voxel, uniform에서 8 voxel 이상입니다. Modal은 voxel 간격을 한 번 더 절반으로 줄였을 때 부피 변화가 3번 `0.071%`, 4번 `0.138%`였습니다. 주기적인 uniform 격자는 5개 sub-voxel offset으로 aliasing을 확인한 뒤 반지름을 각각 `1.015`, `1.005`배 보정했습니다. 해석 시 유효 uniform 반지름은 3번 `0.818196 mm`, 4번 `0.413683 mm`입니다. 이는 부피 수렴 검증이며 주파수 mesh 수렴 검증은 아닙니다.

직접 voxel union의 폐쇄 외벽 STL을 추출하려면:

```powershell
python .\export_voxel_surface.py modal_cantilever4
```

CSV에서 4번 3MF를 다시 만들려면:

```powershell
python .\csv_to_thick_graph_3mf.py
```

## 3D 프린팅용 경량 OBJ

Release의 `*_print.obj.zip`은 원본 대비 face 수를 약 90% 줄이고 MeshFix로 폐쇄·manifold 상태를 복구한 파일입니다.

| 형상 | 원본 face | 경량 face | 부피 변화 |
|---|---:|---:|---:|
| modal | 49,907,540 | 4,987,952 | -0.177% |
| uniform | 48,789,056 | 4,875,214 | -0.184% |

두 경량 모델 사이의 부피 차이는 약 `0.019%`입니다. 다시 만들려면 `python .\reduce_obj_for_print.py`를 실행합니다.
