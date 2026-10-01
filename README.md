# cantilever-new

`modal3`와 `uniform3` 형상을 pygalmesh volume mesh → ANSYS SOLID187 TET10으로 변환하고, 세 재료의 `x < 0` 고정 modal 해석을 실행합니다.

## 입력 형상

OBJ는 용량 때문에 Git 저장소가 아닌 `v1.0` Release에 압축해 두었습니다. `run.py`는 OBJ가 없으면 GitHub CLI로 다음 파일을 자동 다운로드·압축 해제합니다.

- `modal_cantilever3.obj.zip`
- `uniform_cantilever3.obj.zip`

공개 저장소이므로 GitHub 인증은 필요하지 않습니다. 최종 OBJ 실체 부피는 modal `83,319.77 mm³`, uniform `83,341.32 mm³`로 차이는 `0.026%`입니다. uniform 공통 반지름은 `0.806104 mm`입니다.

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

1. OBJ를 pygalmesh TET4 volume mesh로 변환
2. edge 중간 절점을 추가해 SOLID187 TET10 CDB 생성
3. `runs/`에서 Inconel 718, Structural Steel, Formlabs Nylon 12 해석
4. `runs/fixed_xlt0_modal_frequencies_all.csv` 저장

경계조건은 `x < 0`인 모든 노드의 `UX/UY/UZ=0`이며, 탄성 모드 `f1–f6`를 추출합니다.

## 3D 프린팅용 경량 OBJ

Release의 `*_print.obj.zip`은 원본 대비 face 수를 약 90% 줄이고 MeshFix로 폐쇄·manifold 상태를 복구한 파일입니다.

| 형상 | 원본 face | 경량 face | 부피 변화 |
|---|---:|---:|---:|
| modal | 49,907,540 | 4,987,952 | -0.177% |
| uniform | 48,789,056 | 4,875,214 | -0.184% |

두 경량 모델 사이의 부피 차이는 약 `0.019%`입니다. 다시 만들려면 `python .\reduce_obj_for_print.py`를 실행합니다.
